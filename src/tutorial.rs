extern crate piston;
extern crate graphics;
extern crate glutin_window;
extern crate opengl_graphics;

use piston::window::{WindowSettings};
use piston::event_loop::*;
use piston::input::*;
use glutin_window::GlutinWindow as Window;
use opengl_graphics::{ GlGraphics, OpenGL };
use fps_counter::FPSCounter;



pub struct MouseXY {
    x: f64,
    y: f64,
}

pub struct Square {
    x: f64,
    y: f64,
    angle: f64,
    size: f64,
    velocity: (f64, f64),
}

pub struct Circle {
    x: f64,
    y: f64,
    r: f64,
    v: (f64, f64),
}

pub struct WindowSize {
    h: u32,
    w: u32,
}

pub struct App {
    gl: GlGraphics, // OpenGL drawing backend.
    mouse_xy: MouseXY,
    circle: Circle,
    window: WindowSize,
    fps: FPSCounter,
}

const GREEN: [f32; 4] = [0.0, 1.0, 0.0, 1.0];
const RED:   [f32; 4] = [1.0, 0.0, 0.0, 1.0];
const YELLOW:   [f32; 4] = [1.0, 1.0, 0.0, 1.0];
const BLUE:   [f32; 4] = [0.0, 0.0, 1.0, 1.0];

impl App {
    fn render(&mut self, args: &RenderArgs) {
        use graphics::*;
        self.window.h = args.draw_height;
        self.window.w = args.draw_width;

        let earth_bounding_rect = rectangle::square(0.0, 0.0, self.circle.r);
        let sun_bounding_rect = rectangle::square(0.0, 0.0, self.circle.r);

        let Circle { x, y, r, v } = self.circle;
        let earth_x = x * self.window.w as f64;
        let earth_y = y * self.window.h as f64;

        let sun_x = self.window.w as f64 * 0.5;
        let sun_y = self.window.h as f64 * 0.5;

        self.gl.draw(args.viewport(), |c, gl| {
            clear(GREEN, gl);
            let earth_transform = c.transform
                .trans(earth_x, earth_y)
                .trans(-r * 0.5, -r * 0.5);
            let sun_transform = c.transform
                .trans(sun_x, sun_y)
                .trans(-r * 0.5, -r * 0.5);
            ellipse(BLUE, earth_bounding_rect, earth_transform, gl);
            ellipse(YELLOW, sun_bounding_rect, sun_transform, gl);
        })
    }

    fn update(&mut self, args: &UpdateArgs) {
//        self.square.angle = self.square.angle + 4.0 * args.dt;
        // every update, move the square's x/y by its velocity and set the velocity equal to itself plus the acceleration due to gravity * dt
        let Circle { x, y, r, v: (vx, vy) } = self.circle;
        // move the square based on its velocity
        self.circle.x += vx;
        self.circle.y += vy;
        // update the velocity based on gravity
        let d = distance(x, y, 0.5, 0.5);
        let dx = 0.5 - x;
        let dy = 0.5 - y;
        let g = 0.0000001 / d.powi(2);
        let gx =  dx.signum() * g * (dx.abs() / (dx.abs() + dy.abs()));
        let gy =  dy.signum() * g * (dy.abs() / (dx.abs() + dy.abs()));
        self.circle.v.0 += gx;
        self.circle.v.1 += gy;
    }
}

fn clip(min: f64, v: f64, max: f64) -> f64 {
    if min > v {
        min
    } else if max < v {
        max
    } else {
        v
    }
}

fn distance(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    ((x1 - x2).powi(2) + (y1 - y2).powi(2)).sqrt()
}

fn main() {
    let opengl = OpenGL::V3_2;

    let mut window: Window = WindowSettings::new("spinning-circle", [800,800])
//        .fullscreen(true)
        .opengl(opengl)
        .exit_on_esc(true)
        .build()
        .expect("failed to create window");

    let mut app = App {
        gl: GlGraphics::new(opengl),
        mouse_xy: MouseXY { x: 0.0, y: 0.0 },
        circle: Circle { x: 0.1, y: 0.5, r: 50.0, v: (0.00005, 0.00005) },
        window: WindowSize { w: 800, h: 800 },
        fps: FPSCounter::new(),
    };

    let mut event_settings = EventSettings::new();
    event_settings.max_fps = 140;
    event_settings.ups = 240;
    let mut events = Events::new(event_settings);

    while let Some(e) = events.next(&mut window) {
        if let Some(Button::Mouse(button)) = e.press_args() {
            if let MouseButton::Left = button {
                app.circle.x = app.mouse_xy.x as f64 / app.window.w as f64;
                app.circle.y = app.mouse_xy.y as f64 / app.window.h as f64;
                app.circle.v = (0.0005, 0.0005);
            }
        }

        if let Some(Button::Keyboard(Key::A)) = e.press_args() {
            println!("accelerating");
            app.circle.v.0 *= 1.1;
            app.circle.v.1 *= 1.1;
        }

        if let Some(Button::Keyboard(Key::D)) = e.press_args() {
            println!("decelerating");
            app.circle.v.0 *= 0.9;
            app.circle.v.1 *= 0.9;
        }

        e.mouse_cursor(|x, y| {
            app.mouse_xy = MouseXY { x, y };
        });

        if let Some(r) = e.render_args() {
            app.render(&r);
        }

        if let Some(u) = e.update_args() {
            app.update(&u);
        }
    }
}
