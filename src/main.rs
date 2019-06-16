extern crate piston;
extern crate graphics;
extern crate glutin_window;
extern crate opengl_graphics;

use std::f64::consts::PI;
use std::ops::Add;

use find_folder;
use fps_counter::FPSCounter;
use graphics::Context;
use opengl_graphics::{ GlGraphics, OpenGL };
use rand::Rng;
use rayon::prelude::*;
use piston::input::*;
use piston::event_loop::*;
use piston::window::WindowSettings;
use glutin_window::GlutinWindow;
use opengl_graphics::GlyphCache;
use opengl_graphics::TextureSettings;
use piston::window::{AdvancedWindow, Size, Window};


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
}

#[derive(Debug)]
struct Satellite {
    id: usize,
    pos: Position,
    vel: Velocity,
}

impl Satellite {
    fn new(pos: Position, vel: Velocity) -> Satellite {
        Satellite { id: 0, pos, vel }
    }

    fn render(&self, context: &Context, gl: &mut GlGraphics) {
        use graphics::*;
        ellipse(RED, [self.pos.x - 2.0, self.pos.y - 2.0, 4.0, 4.0], context.transform, gl);
    }

    fn update(&self, a: Acceleration) -> Satellite {
        Satellite {
            id: self.id,
            pos: self.vel.apply(self.pos),
            vel: a.apply(self.vel),
        }
    }
    
    fn is_onscreen(&self, w: f64, h: f64) -> bool {
        -2.0 <= self.pos.x && self.pos.x < w + 2.0 && 
        -2.0 <= self.pos.y && self.pos.y < h + 2.0
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

struct State {
    next_planet_id: usize,
    planets: Vec<Planet>,
    next_satellite_id: usize,
    satellites: Vec<Satellite>,
    simulation_speed: usize,
    pause: bool,
    fps: FPSCounter,
    // id of the focused planet
    focused: Option<usize>,
    // id of the moused over planet, used in cases where focused is none
    hovered: Option<usize>,
    cursor: [f64; 2],
}

impl State {
    fn new() -> State {
        State {
            next_planet_id: 0,
            planets: vec![],
            next_satellite_id: 0,
            satellites: vec![],
            simulation_speed: 1,
            pause: false,
            fps: FPSCounter::new(),
            focused: None,
            hovered: None,
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
                Velocity::new(rng.gen_range(-2.0, 2.0), rng.gen_range(-2.0, 2.0)),
                rng.gen_range(150_000_000_000.0, 300_000_000_000.0),
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
                if let Some(focused_pid) = self.focused {
                    for (absorber, absorbed) in colliding {
                        if focused_pid == absorber.id {
                            // focus can stay the same
                            break
                        } else if absorbed.iter().any(|p2| p2.id == focused_pid) {
                            // shift focus to absorber
                            self.focused = Some(absorber.id);
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
                    // check for satellites colliding with planets
                    .filter(|s| {
                        !updated_planets.par_iter().any(|p| s.pos.distance(&p.pos) < p.r)
                    })
                    .collect();

                // apply updates at the end
                self.planets = updated_planets;
                self.satellites = updated_satellites;
            }
        }
    }

    fn render(&mut self, args: &RenderArgs, font: &mut GlyphCache, gl: &mut GlGraphics) {
        // check if the hovered planet has changed
        self.hovered = self.planets.iter().filter(|p| p.contains_point(self.cursor[0], self.cursor[1])).next().map(|p| p.id);

        use graphics::*;
        gl.draw(args.viewport(), |c, gl| {
            clear(BLACK, gl);
            self.planets.iter().for_each(|p| p.render(&c, gl));
            self.satellites.iter().for_each(|s| s.render(&c, gl));

            if let Some(pid) = self.focused.or(self.hovered) {
                let p = self.planets.iter().find(|p| p.id == pid).unwrap();
                let r = p.r * 1.5;
                let ll = p.r * 3.0 / 2.5;
                let x = p.pos.x;
                let y = p.pos.y;
                line(GREEN, 1.0, [x - r, y - r, x - r + ll, y - r], c.transform, gl);
                line(GREEN, 1.0, [x - r, y - r, x - r, y - r + ll], c.transform, gl);
                line(GREEN, 1.0, [x + r, y - r, x + r - ll, y - r], c.transform, gl);
                line(GREEN, 1.0, [x + r, y - r, x + r, y - r + ll], c.transform, gl);
                line(GREEN, 1.0, [x - r, y + r, x - r + ll, y + r], c.transform, gl);
                line(GREEN, 1.0, [x - r, y + r, x - r, y + r - ll], c.transform, gl);
                line(GREEN, 1.0, [x + r, y + r, x + r - ll, y + r], c.transform, gl);
                line(GREEN, 1.0, [x + r, y + r, x + r, y + r - ll], c.transform, gl);
                text(WHITE, 22, &format!("radius: {:.1?}", p.r), font, c.transform.trans(100.0, 300.0), gl).unwrap();
                text(WHITE, 22, &format!("mass: {:.1?}", p.mass), font, c.transform.trans(100.0, 340.0), gl).unwrap();
                text(WHITE, 22, &format!("speed: {:.1?}", p.vel.speed()), font, c.transform.trans(100.0, 380.0), gl).unwrap();
                text(WHITE, 22, &format!("position: {:.1?}, {:.1?}", p.pos.x, p.pos.y), font, c.transform.trans(100.0, 420.0), gl).unwrap();
            }

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
            self.random_setup(500, 0, w, h);
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
                    state.focused = state.hovered;
                },
                Button::Mouse(MouseButton::Right) => {
                    state.planets.clear();
                    state.satellites.clear();
                    let Size { width, height } = window.size();
                    state.random_setup(500, 0, width, height);
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

        e.render(|r| state.render(r, &mut font, &mut gl));
        e.update(|u| state.update(u));
    }
}