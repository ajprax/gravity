extern crate piston;
extern crate graphics;
extern crate glutin_window;
extern crate opengl_graphics;

use std::f64::consts::PI;
use std::ops::Add;

use find_folder;
use fps_counter::FPSCounter;
use glutin_window::GlutinWindow;
use graphics::Context;
use graphics::types::Color;
use graphics::color
use num_format::{Locale, ToFormattedString};
use opengl_graphics::{GlGraphics, GlyphCache, OpenGL, TextureSettings};
use piston::event_loop::*;
use piston::input::*;
use piston::window::{AdvancedWindow, Size, Window, WindowSettings};
use rand::Rng;
use rayon::prelude::*;


/// take a rectangle described by any two points and find the top left (min) and bottom right (max)
/// corners
fn rect_as_min_max(p1: [f64; 2], p2: [f64; 2]) -> [[f64;2]; 2] {
    let [x1, y1] = p1;
    let [x2, y2] = p2;
    let (minx, maxx) = if x1 < x2 { (x1, x2) } else { (x2, x1) };
    let (miny, maxy) = if y1 < y2 { (y1, y2) } else { (y2, y1) };
    [[minx, miny], [maxx, maxy]]
}


fn render_bounding_box(color: Color, min: [f64; 2], max: [f64; 2], c: Context, gl: &mut GlGraphics) {
    use graphics::line;
    let [minx, miny] = min;
    let [maxx, maxy] = max;
    line(color, 1.0, [minx, miny, maxx, miny], c.transform, gl);
    line(color, 1.0, [maxx, miny, maxx, maxy], c.transform, gl);
    line(color, 1.0, [maxx, maxy, minx, maxy], c.transform, gl);
    line(color, 1.0, [minx, maxy, minx, miny], c.transform, gl);
}


#[derive(Clone, Copy, Debug)]
struct Position {
    x: f64,
    y: f64,
}

impl Position {
    fn new(x: f64, y: f64) -> Position {
        Position { x, y }
    }

    fn distance(&self, other: &Position) -> f64 {
        ((self.x - other.x).powi(2) + (self.y - other.y).powi(2)).sqrt()
    }

    fn in_rect(&self, (start, end): ([f64; 2], [f64; 2])) -> bool {
        let [[minx, miny], [maxx, maxy]] = rect_as_min_max(start, end);
        minx <= self.x && self.x < maxx && miny <= self.y && self.y < maxy
    }
}

#[derive(Clone, Copy, Debug)]
struct Velocity {
    x: f64,
    y: f64,
}

impl Velocity {
    fn new(x: f64, y: f64) -> Velocity {
        Velocity { x, y }
    }

    fn apply(&self, pos: Position) -> Position {
        Position::new(self.x + pos.x, self.y + pos.y)
    }

    fn speed(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2)).sqrt()
    }
}

impl Add for Velocity {
    type Output = Velocity;

    fn add(self, rhs: Velocity) -> Velocity {
        Velocity::new(self.x + rhs.x, self.y + rhs.y)
    }
}

#[derive(Clone, Copy, Debug)]
struct Acceleration {
    x: f64,
    y: f64,
}

impl Acceleration {
    fn new(x: f64, y: f64) -> Acceleration {
        Acceleration { x, y }
    }

    fn zero() -> Acceleration {
        Acceleration { x: 0.0, y: 0.0 }
    }

    fn apply(&self, vel: Velocity) -> Velocity {
        Velocity::new(self.x + vel.x, self.y + vel.y)
    }
}

impl Add for Acceleration {
    type Output = Acceleration;

    fn add(self, rhs: Acceleration) -> Acceleration {
        Acceleration::new(self.x + rhs.x, self.y + rhs.y)
    }
}

#[derive(Clone, Debug)]
struct Planet {
    id: usize,
    r: f64,
    pos: Position,
    vel: Velocity,
    mass: f64,
    color: [f32; 4],
}

impl Planet {
    fn new(r: f64, pos: Position, vel: Velocity, mass: f64, color: [f32; 4]) -> Planet {
        Planet { id: 0, r, pos, vel, mass, color }
    }

    fn area(&self) -> f64 {
        PI * self.r * self.r
    }

    fn gravity(&self, pos: &Position) -> Acceleration {
        let d = self.pos.distance(&pos);
        let dx = self.pos.x - pos.x;
        let dy = self.pos.y - pos.y;
        let g = G * self.mass / d.powi(2);
        Acceleration::new(dx.signum() * g * (dx.abs() / (dx.abs() + dy.abs())), dy.signum() * g * (dy.abs() / (dx.abs() + dy.abs())))
    }

    fn absorb(&self, other: &Planet) -> Planet {
        let m = self.mass + other.mass;
        Planet {
            id: self.id,
            r: ((self.area() + other.area()) / PI).sqrt(),
            pos: {
                Position::new(
                    self.pos.x * self.mass / m + other.pos.x * other.mass / m,
                    self.pos.y * self.mass / m + other.pos.y * other.mass / m,
                )
            },
            vel: {
                Velocity::new(
                    self.vel.x * self.mass / m + other.vel.x * other.mass / m,
                    self.vel.y * self.mass / m + other.vel.y * other.mass / m,
                )
            },
            mass: m,
            color: {
                let sm = self.mass as f32;
                let om = other.mass as f32;
                let m = m as f32;
                [
                    self.color[0] * sm / m + other.color[0] * om / m,
                    self.color[1] * sm / m + other.color[1] * om / m,
                    self.color[2] * sm / m + other.color[2] * om / m,
                    self.color[3] * sm / m + other.color[3] * om / m,
                ]
            },
        }
    }

    fn render(&self, c: &Context, gl: &mut GlGraphics) {
        use graphics::*;
        ellipse(self.color, [self.pos.x - self.r, self.pos.y - self.r, self.r * 2.0, self.r * 2.0], c.transform, gl);
    }

    fn update(&self, a: Acceleration) -> Planet {
        Planet {
            id: self.id,
            r: self.r,
            pos: self.vel.apply(self.pos),
            vel: a.apply(self.vel),
            mass: self.mass,
            color: self.color,
        }
    }

    fn is_onscreen(&self, w: f64, h: f64) -> bool {
        -self.r <= self.pos.x && self.pos.x < w + self.r &&
        -self.r <= self.pos.y && self.pos.y < h + self.r
    }

    fn contains_point(&self, x: f64, y: f64) -> bool {
        self.pos.distance(&Position::new(x, y)) < self.r
    }

    fn facts(&self) -> Vec<String> {
        vec![
            format!("mass: {}", (self.mass as usize).to_formatted_string(&Locale::en)),
            format!("radius: {:.2?}", self.r),
            format!("position: {:.2?}, {:.2?}", self.pos.x, self.pos.y),
            format!("velocity: {:.2?}, {:.2?}", self.vel.x, self.vel.y),
        ]
    }
}

#[derive(Debug)]
struct Satellite {
    id: usize,
    pos: Position,
    vel: Velocity,
}

impl Satellite {
    const R: f64 = 2.0;

    fn new(pos: Position, vel: Velocity) -> Satellite {
        Satellite { id: 0, pos, vel }
    }

    fn render(&self, context: &Context, gl: &mut GlGraphics) {
        use graphics::*;
        ellipse(RED, [self.pos.x - Satellite::R, self.pos.y - Satellite::R, Satellite::R * 2.0, Satellite::R * 2.0], context.transform, gl);
    }

    fn update(&self, a: Acceleration) -> Satellite {
        Satellite {
            id: self.id,
            pos: self.vel.apply(self.pos),
            vel: a.apply(self.vel),
        }
    }
    
    fn is_onscreen(&self, w: f64, h: f64) -> bool {
        -Satellite::R <= self.pos.x && self.pos.x < w + Satellite::R &&
        -Satellite::R <= self.pos.y && self.pos.y < h + Satellite::R
    }

    fn contains_point(&self, x: f64, y: f64) -> bool {
        self.pos.distance(&Position::new(x, y)) < Satellite::R
    }

    fn facts(&self) -> Vec<String> {
        vec![
            format!("position: {:.2?}, {:.2?}", self.pos.x, self.pos.y),
            format!("velocity: {:.2?}, {:.2?}", self.vel.x, self.vel.y),
        ]
    }
}

const BLACK:  [f32; 4] = [0.0, 0.0, 0.0, 0.0];
const GREEN:  [f32; 4] = [0.0, 1.0, 0.0, 1.0];
const RED:    [f32; 4] = [1.0, 0.0, 0.0, 1.0];
const YELLOW: [f32; 4] = [1.0, 1.0, 0.0, 1.0];
const BLUE:   [f32; 4] = [0.0, 0.0, 1.0, 1.0];
const WHITE:  [f32; 4] = [1.0, 1.0, 1.0, 1.0];
const COLORS: [[f32; 4]; 4] = [GREEN, RED, YELLOW, BLUE];

fn random_color() -> [f32; 4]{
    let mut rng = rand::thread_rng();
    [
        rng.gen_range(0.0, 1.0),
        rng.gen_range(0.0, 1.0),
        rng.gen_range(0.0, 1.0),
        rng.gen_range(0.0, 1.0)
    ]
}

const G: f64 = 6.67408e-11;

// the id of the focused (or hovered) object
#[derive(Copy, Clone, Debug)]
enum Focus {
    Planet(usize),
    Satellite(usize),
}

struct State {
    next_planet_id: usize,
    planets: Vec<Planet>,
    next_satellite_id: usize,
    satellites: Vec<Satellite>,
    simulation_speed: usize,
    pause: bool,
    fps: FPSCounter,
    focused: Option<Focus>,
    hovered: Option<Focus>,
    drag_select_start: Option<[f64; 2]>,
    cursor: [f64; 2],
}

impl State {
    fn new() -> State {
        State {
            next_planet_id: 0,
            planets: vec![],
            next_satellite_id: 0,
            satellites: vec![],
            simulation_speed: 5,
            pause: false,
            fps: FPSCounter::new(),
            focused: None,
            hovered: None,
            drag_select_start: None,
            cursor: [0.0, 0.0],
        }
    }

    fn add_planet(&mut self, mut planet: Planet) {
        // when State takes ownership of a planet it overwrites whatever ID was there
        let id = self.next_planet_id;
        self.next_planet_id += 1;
        planet.id = id;
        self.planets.push(planet);
    }

    fn add_satellite(&mut self, mut satellite: Satellite) {
        // when State takes ownership of a satellite it overwrites whatever ID was there
        let id = self.next_satellite_id;
        self.next_satellite_id += 1;
        satellite.id = id;
        self.satellites.push(satellite);
    }

    fn simple_orbit_setup(&mut self) {
        self.add_planet(Planet::new(
            25.0,
            Position::new(1280.0, 220.0),
            Velocity::new(0.0, 0.0),
            500_000_000_000.0,
            BLUE,
        ));

        self.add_satellite(
            Satellite::new(
                Position::new(1330.0, 550.0),
                Velocity::new(0.3, 0.0),
            )
        )
    }

    fn orbital_velocity_experiment_setup(&mut self) {
        self.add_planet(Planet::new(
            100.0,
            Position::new(1280.0, 720.0),
            Velocity::new(0.0, 0.0),
            500_000_000_000.0,
            BLUE,
        ));

        for v in 200..700 {
            self.add_satellite(Satellite::new(
                Position::new(1280.0, 550.0),
                Velocity::new(v as f64 / 1000.0, 0.0),
            ));
        }
    }

    fn rogue_planet_setup(&mut self) {
        self.orbital_velocity_experiment_setup();

        self.add_planet(Planet::new(
            25.0,
            Position::new(1000.0, 550.0),
            Velocity::new(0.6, 0.0),
            1_000_000_000_000.0,
            BLUE,
        ));
        self.add_planet(Planet::new(
            25.0,
            Position::new(1500.0, 950.0),
            Velocity::new(0.1, -0.8),
            1_000_000_000_000.0,
            BLUE,
        ));
    }

    fn random_setup(&mut self, num_planets: usize, num_satellites: usize, w: f64, h: f64) {
        let mut rng = rand::thread_rng();
        for _ in 0..num_planets {
            self.add_planet(Planet::new(
                rng.gen_range(5.0, 10.0),
                Position::new(rng.gen_range(0.0, w), rng.gen_range(0.0, h)),
                Velocity::new(rng.gen_range(-0.25, 0.25), rng.gen_range(-0.25, 0.25)),
                rng.gen_range(1_500_000_000.0, 3_000_000_000.0),
                random_color(),
            ));
        }
        for _ in 0..num_satellites {
            self.add_satellite(
                Satellite::new(
                    Position::new(rng.gen_range(0.0, w), rng.gen_range(0.0, h)),
                    Velocity::new(rng.gen_range(-1.0, 1.0), rng.gen_range(-1.0, 1.0)),
                )
            )
        }
    }

    fn increase_simulation_speed(&mut self) {
        self.simulation_speed += 1;
    }

    fn decrease_simulation_speed(&mut self) {
        if self.simulation_speed > 1 {
            self.simulation_speed -= 1;
        }
    }

    fn update(&mut self, args: &UpdateArgs) {
        // be sure to calculate all updates before applying any of them
        if !self.pause {
            // faster simulation speed means we run the update multiple times so that we retain the precision from smaller steps
            for _ in 0..self.simulation_speed {
                // for each planet and satellite, calculate the effect of gravity of all (other) planets,
                // add them together, and apply them
                let updated_planets: Vec<Planet> = self.planets.par_iter().map(|p| {
                    let a = self.planets.par_iter()
                        .filter(|p2| p.id != p2.id)
                        .map(|p2| p2.gravity(&p.pos))
                        .reduce(|| Acceleration::zero(), |a1, a2| a1 + a2);
                    p.update(a)
                }).collect();
                // check for planetary collisions, detect collisions first then resolve them
                // TODO: handle case of A colliding with B colliding with C, this is pretty unlikely to happen within one update
                let colliding: Vec<(Planet, Vec<&Planet>)> =  updated_planets.par_iter().enumerate().filter_map(|(i, p)| {
                    let colliding: Vec<(usize, &Planet)> = updated_planets.iter().filter(|p2| p.id != p2.id).enumerate().filter_map(|(i2, p2)| {
                        // planets are considered colliding if the center of one is within another
                        if p.pos.distance(&p2.pos) < p.r.max(p2.r) {
                            Some((i2, p2))
                        } else {
                            None
                        }
                    }).collect();

                    if colliding.len() == 0 {
                        // no collisions
                        Some((p.clone(), colliding.into_iter().map(|(i, p2)| p2).collect()))
                    } else if colliding.par_iter().any(|(i2, p2)| *i2 < i) {
                        // this planet was already absorbed by one earlier in the list, so just drop it
                        None
                    } else{
                        // leftmost participant in a collision, absorb the others
                        Some((p.clone(), colliding.into_iter().map(|(i, p2)| p2).collect()))
                    }
                }).collect();
                // resolve collisions
                let updated_planets: Vec<Planet> = colliding.par_iter().map(|(p, to_absorb)| {
                    to_absorb.into_iter().fold(p.clone(), |p1, p2| p1.absorb(p2))
                }).collect();

                // if the focused planet got absorbed, move focus to the absorbing planet
                if let Some(Focus::Planet(focused_pid)) = self.focused {
                    for (absorber, absorbed) in colliding {
                        if focused_pid == absorber.id {
                            // focus can stay the same
                            break
                        } else if absorbed.iter().any(|p2| p2.id == focused_pid) {
                            // shift focus to absorber
                            self.focused = Some(Focus::Planet(absorber.id));
                            break
                        }
                    }
                }

                let updated_satellites: Vec<Satellite> = self.satellites.par_iter()
                    .map(|s| {
                        let a = self.planets.par_iter()
                            .map(|p| p.gravity(&s.pos))
                            .reduce(Acceleration::zero, |a1, a2| a1 + a2);
                        s.update(a)
                    })
                    .collect();
                // check for satellites colliding with planets
                let (updated_satellites, crashing_satellites): (Vec<Satellite>, Vec<Satellite>) = updated_satellites.into_iter().partition(|s| {
                    !updated_planets.par_iter().any(|p| s.pos.distance(&p.pos) < p.r)
                });
                // if the focused satellite crashes into a planet, unfocus it
                if let Some(Focus::Satellite(sid)) = self.focused {
                    if crashing_satellites.iter().any(|s| s.id == sid) {
                        self.focused = None;
                    }
                }
                // apply updates at the end
                self.planets = updated_planets;
                self.satellites = updated_satellites;
            }
        }
    }

    fn render(&mut self, args: &RenderArgs, font: &mut GlyphCache, gl: &mut GlGraphics) {
        // check if the hovered target has changed, but only if we're not drag selecting
        // checked in render because hovered only used in render so we don't need to pay the cost on
        // every update
        if self.drag_select_start.is_none() {
            self.hovered = self.planets.iter().find(|p| p.contains_point(self.cursor[0], self.cursor[1])).map(|p| Focus::Planet(p.id)).or(
                self.satellites.iter().find(|s| s.contains_point(self.cursor[0], self.cursor[1])).map(|s| Focus::Satellite(s.id))
            );
        }

        use graphics::*;
        gl.draw(args.viewport(), |c, gl| {
            clear(BLACK, gl);
            // draw planets and satellites
            self.planets.iter().for_each(|p| p.render(&c, gl));
            self.satellites.iter().for_each(|s| s.render(&c, gl));
            // draw the box selection
            // TODO: preview which item will be selected?
            if let Some(start) = self.drag_select_start {
                let [min, max] = rect_as_min_max(start, self.cursor);
                render_bounding_box(GREEN, min, max, c, gl);
            }
            // draw details about the selected planet or satellite
            let deets: Option<(f64, f64, f64, f64, Vec<String>)> = match self.focused.or(self.hovered) {
                Some(Focus::Planet(pid)) => {
                    let p = self.planets.iter().find(|p| p.id == pid).unwrap();
                    let r = p.r * 1.5;
                    let ll = p.r * 3.0 / 2.5;
                    let x = p.pos.x;
                    let y = p.pos.y;
                    Some((r, ll, x, y, p.facts()))
                },
                Some(Focus::Satellite(sid)) => {
                    let s = self.satellites.iter().find(|s| s.id == sid).unwrap();
                    let r = Satellite::R * 1.5;
                    let ll = Satellite::R * 3.0 / 2.5;
                    let x = s.pos.x;
                    let y = s.pos.y;
                    Some((r, ll, x, y, s.facts()))
                },
                None => None,
            };
            if let Some((r, ll, x, y, facts)) = deets {
                line(GREEN, 1.0, [x - r, y - r, x - r + ll, y - r], c.transform, gl);
                line(GREEN, 1.0, [x - r, y - r, x - r, y - r + ll], c.transform, gl);
                line(GREEN, 1.0, [x + r, y - r, x + r - ll, y - r], c.transform, gl);
                line(GREEN, 1.0, [x + r, y - r, x + r, y - r + ll], c.transform, gl);
                line(GREEN, 1.0, [x - r, y + r, x - r + ll, y + r], c.transform, gl);
                line(GREEN, 1.0, [x - r, y + r, x - r, y + r - ll], c.transform, gl);
                line(GREEN, 1.0, [x + r, y + r, x + r - ll, y + r], c.transform, gl);
                line(GREEN, 1.0, [x + r, y + r, x + r, y + r - ll], c.transform, gl);
                let mut fact_y = 300.0;
                for fact in facts {
                    text(WHITE, 22, &fact, font, c.transform.trans(100.0, fact_y), gl).unwrap();
                    fact_y += 40.0;
                }
            }
            // draw general details about the state of the world
            text(WHITE, 22, &format!("fps: {}", self.fps.tick()), font, c.transform.trans(100.0, 100.0), gl).unwrap();
            text(WHITE, 22, &format!("planets: {}", self.planets.len()), font, c.transform.trans(100.0, 150.0), gl).unwrap();
            text(WHITE, 22, &format!("speed: {}x{}", self.simulation_speed, if self.pause { " (paused)" } else { "" }), font, c.transform.trans(100.0, 200.0), gl).unwrap();
        });

        let [w, h] = args.viewport().window_size;
        if self.planets.len() == 1 ||
            (!self.planets.par_iter().any(|p| p.is_onscreen(w, h)) &&
             !self.satellites.par_iter().any(|s| s.is_onscreen(w, h)))
        {
            self.planets.clear();
            self.satellites.clear();
            self.focused = None;
            self.random_setup(500, 20, w, h);
        }
    }
}

fn main() {
    let opengl = OpenGL::V4_5;

    let mut window: GlutinWindow = WindowSettings::new("gravity", [0, 0])
        .graphics_api(opengl)
        .fullscreen(true)
        .exit_on_esc(true)
        .build()
        .expect("failed to create window");

    let mut font = GlyphCache::new(
        find_folder::Search::ParentsThenKids(3, 3).for_folder("assets").unwrap().join("Roboto-Regular.ttf"),
        (),
        TextureSettings::new()
    ).unwrap();

    let mut gl = GlGraphics::new(opengl);
    let mut state = State::new();

    let mut event_settings = EventSettings::new();
    event_settings.max_fps = 140;
    event_settings.ups = 240;
    let mut events = Events::new(event_settings);

    while let Some(e) = events.next(&mut window) {
        // we don't know what size the screen is, but with fullscreen(true) we'll get a resize event
        // that sets us to the full size of the screen before getting an event that sets us to the
        // configured size. we take advantage of that event to set the size correctly.
        e.resize(|r| {
            let Size { width, height } = window.size();
            let [w, h] = r.window_size;
            if w > width || h > height {
                window.set_size([w, h]);
            }
        });

        e.mouse_cursor(|c| {
            state.cursor = c;
        });

        e.press(|b| {
            match b {
                Button::Mouse(MouseButton::Left) => {
                    state.drag_select_start = Some(state.cursor);
                    // while drag selecting, nothing should be considered hovered
                    state.hovered = None;
                },
                Button::Mouse(MouseButton::Right) => {
                    state.planets.clear();
                    state.satellites.clear();
                    let Size { width, height } = window.size();
                    state.random_setup(500, 20, width, height);
                },
                Button::Keyboard(Key::NumPadPlus) | Button::Keyboard(Key::Equals) => {
                    state.increase_simulation_speed();
                },
                Button::Keyboard(Key::NumPadMinus) | Button::Keyboard(Key::Minus) => {
                    state.decrease_simulation_speed();
                },
                Button::Keyboard(Key::Space) => {
                    state.pause = !state.pause;
                },
                _ => ()
            }
        });

        e.release(|b| {
            match b {
                Button::Mouse(MouseButton::Left) => {
                    if let Some(start) = state.drag_select_start {
                        state.drag_select_start = None;
                        let end = state.cursor;
                        // TODO: multiselect?
                        for p in &state.planets {
                            if p.pos.in_rect((start, end)) || p.contains_point(end[0], end[1]) {
                                state.focused = Some(Focus::Planet(p.id));
                                return
                            }
                        }
                        for s in &state.satellites {
                            if s.pos.in_rect((start, end)) || s.contains_point(end[0], end[1]) {
                                state.focused = Some(Focus::Satellite(s.id));
                                return
                            }
                        }
                        state.focused = None;
                    }
                },
                _ => ()
            }
        });

        e.render(|r| state.render(r, &mut font, &mut gl));
        e.update(|u| state.update(u));
    }
}