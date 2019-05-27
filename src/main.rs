extern crate piston;
extern crate graphics;
extern crate glutin_window;
extern crate opengl_graphics;

use std::f64::consts::PI;
use std::ops::Add;

use fps_counter::FPSCounter;
use graphics::Context;
use piston::window::WindowSettings;
use piston::event_loop::*;
use piston::input::*;
use glutin_window::GlutinWindow as Window;
use opengl_graphics::{ GlGraphics, OpenGL };



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

    fn apply(&self, loc: Position) -> Position {
        Position::new(self.x + loc.x, self.y + loc.y)
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
}

impl Planet {
    fn new(r: f64, pos: Position, vel: Velocity, mass: f64) -> Planet {
        Planet { id: 0, r, pos, vel, mass }
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
        }
    }

    fn render(&self, context: &Context, gl: &mut GlGraphics) {
        use graphics::*;
        ellipse(BLUE, [self.pos.x - self.r, self.pos.y - self.r, self.r * 2.0, self.r * 2.0], context.transform, gl);
    }

    fn update(&self, a: Acceleration) -> Planet {
        Planet {
            id: self.id,
            r: self.r,
            pos: self.vel.apply(self.pos),
            vel: a.apply(self.vel),
            mass: self.mass,
        }
    }
}

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
}

const BLACK:  [f32; 4] = [0.0, 0.0, 0.0, 0.0];
const GREEN:  [f32; 4] = [0.0, 1.0, 0.0, 1.0];
const RED:    [f32; 4] = [1.0, 0.0, 0.0, 1.0];
const YELLOW: [f32; 4] = [1.0, 1.0, 0.0, 1.0];
const BLUE:   [f32; 4] = [0.0, 0.0, 1.0, 1.0];
const G: f64 = 6.67408e-11;

struct App {
    gl: GlGraphics,
    fps: FPSCounter,
}

impl App {
    fn new(gl: GlGraphics) -> App {
        App { gl, fps: FPSCounter::new() }
    }

    fn render(&mut self, args: &RenderArgs, state: &State) {
        use graphics::clear;
        self.gl.draw(args.viewport(), |c, gl| {
            clear(BLACK, gl);
            state.planets.iter().for_each(|p| p.render(&c, gl));
            state.satellites.iter().for_each(|s| s.render(&c, gl));
        });
    }

    fn update(&mut self, args: &UpdateArgs, state: &mut State) {
        // be sure to calculate all updates before applying any of them

        // for each planet and satellite, calculate the effect of gravity of all (other) planets,
        // add them together, and apply them
        let updated_planets: Vec<Planet> = state.planets.iter().map(|p| {
            let a = state.planets.iter()
                .filter(|p2| p.id != p2.id)
                .map(|p2| p2.gravity(&p.pos))
                .fold(Acceleration::zero(), |a1, a2| a1 + a2);
            p.update(a)
        }).collect();
        // check for planetary collisions
        let updated_planets: Vec<Planet> = updated_planets.iter().enumerate().filter_map(|(i, p)| {
            let colliding: Vec<(usize, &Planet)> = updated_planets.iter().filter(|p2| p.id != p2.id).enumerate().filter_map(|(i2, p2)| {
                // planets are considered colliding if the center of one is within another
                if p.pos.distance(&p2.pos) < (p.r - p2.r).abs() {
                    Some((i2, p2))
                } else {
                    None
                }
            }).collect();

            if colliding.len() == 0 {
                // no collisions
                Some(p.clone())
            } else if colliding.iter().any(|(i2, p2)| i2 < &i) {
                // this planet was already absorbed by one earlier in the list, so just drop it
                None
            } else{
                // leftmost participant in a collision, absorb the others
                let absorbed: Planet = colliding.into_iter().fold(p.clone(), |p1, (i2, p2)| p1.absorb(p2));
                Some(absorbed)
            }
            // TODO: handle case of A colliding with B colliding with C, this is pretty unlikely to happen within one update
        }).collect();
        let updated_satellites: Vec<Satellite> = state.satellites.iter()
            .map(|s| {
                let a = state.planets.iter()
                    .map(|p| p.gravity(&s.pos))
                    .fold(Acceleration::zero(), |a1, a2| a1 + a2);
                s.update(a)
            })
            // check for satellites colliding with planets
            .filter(|s| {
                !updated_planets.iter().any(|p| s.pos.distance(&p.pos) < p.r)
            })
            .collect();
        // apply updates at the end
        state.planets = updated_planets;
        state.satellites = updated_satellites;
    }
}

struct State {
    next_planet_id: usize,
    planets: Vec<Planet>,
    next_satellite_id: usize,
    satellites: Vec<Satellite>,
}

impl State {
    fn new() -> State {
        State { next_planet_id: 0, planets: vec![], next_satellite_id: 0, satellites: vec![] }
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
}

fn orbital_velocity_experiment_setup(state: &mut State) {
    state.add_planet(Planet::new(
        100.0,
        Position::new(1280.0, 720.0),
        Velocity::new(0.0, 0.0),
        500_000_000_000.0
    ));

    for v in 200..700 {
        state.add_satellite(Satellite::new(
            Position::new(1280.0, 550.0),
            Velocity::new(v as f64 / 1000.0, 0.0),
        ));
    }
}


fn rogue_planet_setup(state: &mut State) {
    orbital_velocity_experiment_setup(state);

    state.add_planet(Planet::new(
        25.0,
        Position::new(1000.0, 550.0),
        Velocity::new(0.6, 0.0),
        1_000_000_000_000.0,
    ));
    state.add_planet(Planet::new(
        25.0,
        Position::new(1500.0, 950.0),
        Velocity::new(0.1, -0.8),
        1_000_000_000_000.0,
    ));
}

fn main() {
    let opengl = OpenGL::V3_2;

    let mut window: Window = WindowSettings::new("gravity", [1000, 1000])
        .opengl(opengl)
        .fullscreen(true)
        .exit_on_esc(true)
        .build()
        .expect("failed to create window");

    let mut app = App::new(GlGraphics::new(opengl));
    let mut state = State::new();

    rogue_planet_setup(&mut state);

    let mut event_settings = EventSettings::new();
    event_settings.max_fps = 140;
    event_settings.ups = 240;
    let mut events = Events::new(event_settings);

    while let Some(e) = events.next(&mut window) {
        if let Some(r) = e.render_args() {
            app.render(&r, &state);
        }

        if let Some(u) = e.update_args() {
            app.update(&u, &mut state);
        }
    }
}