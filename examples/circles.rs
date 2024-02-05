#![allow(dead_code)]
use std::collections::{HashMap, VecDeque};

use circles_lines::set_center;
use macroquad::prelude::*;

#[derive(Debug)]
pub struct Circle {
    pub pos: Vec2,
    pub vel: Vec2,
    pub acc: Vec2,
    pub radius: f32,
}

impl Circle {
    pub fn new(pos: Vec2, vel: Vec2) -> Self {
        Self {
            pos,
            vel,
            acc: Vec2::ZERO,
            radius: 2.0,
        }
    }

    pub fn update(&mut self, dt: f32) {
        self.acc += vec2(0.0, 0.98); // Gravity
        self.vel += std::mem::take(&mut self.acc) * dt;
        self.vel *= 0.999f32.powf(dt * 60.0f32); // Every frame we mul by 0.999
        self.pos += self.vel * dt;
    }

    pub fn draw(&self) {
        draw_circle(self.pos.x, self.pos.y, self.radius, WHITE);
    }

    pub fn draw_vel(&self) {
        let (x1, y1) = self.pos.into();
        let (x2, y2) = (self.pos + self.vel * 2.0).into();
        draw_line(x1, y1, x2, y2, 2.0, RED);
    }
}

#[derive(Debug)]
pub struct Edge {
    pub start: usize,
    pub end: usize,
    pub rest_distance: f32,
}

impl Edge {
    pub fn new(start: usize, end: usize, rest_distance: f32) -> Self {
        Self {
            start,
            end,
            rest_distance,
        }
    }
}

#[derive(Debug)]
pub struct MouseState {
    pub prev_pos: Vec2,
    pub prev_vel: Vec2,
    pub prev_is_down: bool,
    pub pos: Vec2,
    pub vel: Vec2,
    pub is_down: bool,
    pub track: VecDeque<Vec2>,
}

impl MouseState {
    pub fn new() -> Self {
        let mut self_ = Self {
            prev_pos: Vec2::ZERO,
            prev_vel: Vec2::ZERO,
            prev_is_down: false,
            pos: Vec2::ZERO,
            vel: Vec2::ZERO,
            is_down: false,
            track: VecDeque::new(),
        };
        self_.handle();
        self_
    }

    pub fn handle(&mut self) {
        let pos: Vec2 = mouse_position().into();
        let vel: Vec2 = mouse_delta_position().into();
        let is_down = is_mouse_button_down(MouseButton::Left);
        {
            self.prev_pos = self.pos;
            self.prev_vel = self.vel;
            self.prev_is_down = self.is_down;
            self.pos = pos;
            self.vel = vel;
            self.is_down = is_down;

            self.track.push_front(pos);
            self.track.truncate(10);
        }
    }

    pub fn is_click(&self) -> bool {
        !self.prev_is_down && self.is_down
    }

    pub fn is_release(&self) -> bool {
        self.prev_is_down && !self.is_down
    }
}

#[derive(Debug)]
pub struct BoudaryState {
    pub xl: f32,
    pub xr: f32,
    pub yt: f32,
    pub yb: f32,
}

impl BoudaryState {
    pub fn new() -> Self {
        let mut self_ = Self {
            xl: 0.0,
            xr: 1280.0,
            yt: 0.0,
            yb: 720.0,
        };
        self_.handle();
        self_
    }
    pub fn handle(&mut self) {
        let w = screen_width();
        let h = screen_height();
        *self = Self {
            xl: 0.0,
            xr: w,
            yt: 0.0,
            yb: h,
        };
    }
}

#[derive(Debug)]
pub struct State {
    pub circles: Vec<Circle>,
    pub edges: Vec<Edge>,
    pub constraints: Vec<(usize, Vec2)>,
    pub touched_circle: Option<usize>,
    pub mouse: MouseState,
    pub boudary: BoudaryState,
    pub dt: f32,
    pub n_updates: usize,
}

pub fn get_line(x: usize, dx: f32) -> (Vec<Circle>, Vec<Edge>) {
    let mut circles = Vec::with_capacity(x);
    let mut edges = Vec::new();

    for ix in 0..x {
        circles.push(Circle::new(
            vec2(ix as f32 * dx + 100.0, 200.0),
            vec2(0.0, 0.0),
        ))
    }
    let distance = |a: usize, b: usize| circles[a].pos.distance(circles[b].pos) / 10.0;
    let edge = |a: usize, b: usize| Edge::new(a, b, distance(a, b));
    for ix in 1..x {
        edges.push(edge(ix - 1, ix));
    }

    (circles, edges)
}

pub fn get_grid(x: usize, y: usize, dx: f32) -> (Vec<Circle>, Vec<Edge>) {
    let mut circles = Vec::with_capacity(x * y);
    let mut edges = Vec::new();
    let mut coords = HashMap::new();
    for ix in 0..x {
        for iy in 0..y {
            coords.insert((ix, iy), circles.len());
            circles.push(Circle::new(
                vec2(ix as f32 * dx + 100.0, iy as f32 * dx + 100.0),
                vec2(0.0, 0.0),
            ));
        }
    }
    let distance = |a: usize, b: usize| circles[a].pos.distance(circles[b].pos);
    let edge = |a: usize, b: usize| Edge::new(a, b, distance(a, b));
    for ix in 0..x {
        for iy in 0..y {
            let i = coords[&(ix, iy)];
            if let Some(&r) = coords.get(&(ix + 1, iy)) {
                edges.push(edge(i, r));
            }
            if let Some(&b) = coords.get(&(ix, iy + 1)) {
                edges.push(edge(i, b));
            }
            if let Some(&rb) = coords.get(&(ix + 1, iy + 1)) {
                edges.push(edge(i, rb));
            }
            if let Some(&rt) = iy.checked_sub(1).and_then(|tiy| coords.get(&(ix + 1, tiy))) {
                edges.push(edge(i, rt));
            }
        }
    }
    let lt = coords[&(0, 0)];
    let rt = coords[&(x - 1, 0)];
    let lb = coords[&(0, y - 1)];
    let rb = coords[&(x - 1, y - 1)];

    // diagonals
    edges.push(edge(lt, rb));
    edges.push(edge(lb, rt));

    // sides
    edges.push(edge(lt, rt)); // top
    edges.push(edge(lt, lb)); // left
    edges.push(edge(lb, rb)); // bottom
    edges.push(edge(rt, rb)); // right

    (circles, edges)
}

#[derive(PartialEq)]
pub enum Example {
    Line,
    Grid,
    ThreePoints,
}

impl State {
    pub fn new() -> Self {
        let example = Example::Line;

        let (mut circles, edges) = match example {
            Example::Line => get_line(40, 25.0),
            Example::Grid => get_grid(6, 6, 60.0),
            Example::ThreePoints => {
                let circles = vec![
                    Circle::new(vec2(100.0, 200.0), vec2(200.0, 0.0)),
                    Circle::new(vec2(200.0, 300.0), vec2(-200.0, -100.0)),
                    Circle::new(vec2(300.0, 100.0), vec2(300.0, 200.0)),
                ];
                let distance = |a: usize, b: usize| circles[a].pos.distance(circles[b].pos);
                let edge = |a: usize, b: usize| Edge::new(a, b, distance(a, b));
                let edges = vec![edge(0, 1), edge(1, 2)];
                (circles, edges)
            }
        };

        let mut positions: Vec<_> = circles.iter_mut().map(|c| &mut c.pos).collect();
        set_center(
            positions.as_mut_slice(),
            vec2(screen_width() / 2.0, screen_height() / 4.0 * 1.0),
        );

        let constraints = vec![
            (0, circles[0].pos),
            (circles.len() - 1, circles.last().unwrap().pos),
        ];

        let n_updates = std::env::args()
            .nth(1)
            .and_then(|s| s.parse().ok())
            .unwrap_or(100);
        Self {
            circles,
            edges,
            constraints,
            touched_circle: None,
            mouse: MouseState::new(),
            boudary: BoudaryState::new(),
            dt: get_frame_time(),
            n_updates,
        }
    }
    pub fn handle(&mut self) {
        self.mouse.handle();
        self.boudary.handle();
    }
    pub fn update(&mut self) {
        for _ in 0..self.n_updates {
            // acceleration
            self.edges.iter().for_each(
                |&Edge {
                     start,
                     end,
                     rest_distance,
                 }| {
                    let (c1, c2) = get_2_mut(&mut self.circles, start, end);
                    let p1 = c1.pos;
                    let p2 = c2.pos;

                    let dist = p1.distance(p2);
                    let acc_mod =
                        (rest_distance - dist).powi(2) * (rest_distance - dist).signum() / 100.0; // C2->acc
                    let n12 = p2 - p1; // C1 -> C2
                    c2.acc += n12.normalize() * acc_mod;
                    c1.acc += -n12.normalize() * acc_mod;
                },
            );
            permutations(self.circles.len()).for_each(|(i1, i2)| {
                let (c1, c2) = get_2_mut(&mut self.circles, i1, i2);
                let p1 = c1.pos;
                let p2 = c2.pos;
                let dist = p1.distance(p2);
                if dist < c1.radius + c2.radius {
                    let acc_mod = 100.0;
                    let n12 = p2 - p1; // C1 -> C2
                    c2.acc += n12.normalize() * acc_mod;
                    c1.acc += -n12.normalize() * acc_mod;
                }
            });
            // velocity
            self.circles.iter_mut().for_each(|c| c.update(self.dt));
            if let Some(touched_circle) = self.touched_circle {
                let circle = &mut self.circles[touched_circle];
                let strength = 0.3;
                circle.pos = Vec2::lerp(circle.pos, self.mouse.pos, strength);

                let strength = 0.8;
                let drag_vel =
                    convert_vel_to_screen(-(self.mouse.vel + self.mouse.prev_vel) / 2.0) * 1.0;
                circle.vel = Vec2::lerp(circle.vel, drag_vel, strength);
                if drag_vel.length() < 5.0 {
                    circle.vel *= 0.0001;
                }
            }
            // boundary
            self.circles.iter_mut().for_each(|c| {
                let BoudaryState { xl, xr, yt, yb } = self.boudary;

                let pos = &mut c.pos;
                let vel = &mut c.vel;

                let r = c.radius;
                let cl = pos.x - r;
                let cr = pos.x + r;
                let ct = pos.y - r;
                let cb = pos.y + r;
                if cl < xl {
                    let dist = (xl - cl).abs();
                    pos.x = xl + r + dist;
                    vel.x *= -1.0;
                }
                if cr > xr {
                    let dist = (xr - cr).abs();
                    pos.x = xr - r - dist;
                    vel.x *= -1.0;
                }
                if ct < yt {
                    let dist = (yt - ct).abs();
                    pos.y = yt + r + dist;
                    vel.y *= -1.0;
                }
                if cb > yb {
                    let dist = (yb - cb).abs();
                    pos.y = yb - r - dist;
                    vel.y *= -1.0;
                }
            });
            // constraints
            self.constraints.iter().for_each(|&(i, pos)| {
                self.circles[i].pos = pos;
                self.circles[i].vel = vec2(0.0, 0.0);
            });
        }
    }
    pub fn draw(&self) {
        self.edges.iter().for_each(|Edge { start, end, .. }| {
            let (x1, y1) = self.circles[*start].pos.into();
            let (x2, y2) = self.circles[*end].pos.into();
            draw_line(x1, y1, x2, y2, 3.0, LIGHTGRAY);
        });
        self.circles.iter().for_each(|c| c.draw());
        self.circles.iter().for_each(|c| c.draw_vel());
    }
}

#[macroquad::main(cfg)]
async fn main() {
    let mut state = State::new();
    loop {
        // Handle inputs
        // keys
        #[cfg(not(target_arch = "wasm32"))]
        if handle_exit() {
            break;
        }
        // mouse
        state.handle();
        handle_touch(&mut state);
        // time
        state.dt = get_frame_time();

        // Update state
        state.update();

        // Draw
        state.draw();

        next_frame().await;
    }
}

fn handle_touch(state: &mut State) {
    let mouse = &state.mouse;
    if mouse.is_click() {
        debug_assert!(state.touched_circle.is_none());
        if let Some(touched_circle) = state.circles.iter().enumerate().find_map(|(i, c)| {
            let test = |mouse_pos: Vec2| mouse_pos.distance(c.pos) < c.radius + 20.0;
            (test(mouse.pos) && mouse.track.iter().copied().any(test)).then_some(i)
        }) {
            state.touched_circle = Some(touched_circle);
        }
    } else if mouse.is_release() {
        if let Some(touched_circle) = state.touched_circle {
            let circle = &mut state.circles[touched_circle];
            let drag_vel = convert_vel_to_screen(-(mouse.vel + mouse.prev_vel) / 2.0) * 1.0;
            circle.vel = drag_vel;
        }
        state.touched_circle = None;
    }
}

/// Convert a position in pixels to a position in the range [-1; 1].
fn convert_to_local(pixel_pos: Vec2) -> Vec2 {
    Vec2::new(pixel_pos.x / screen_width(), pixel_pos.y / screen_height()) * 2.0
        - Vec2::new(1.0, 1.0)
}

/// Convert a position in range [-1; 1] to a position in pixels.
fn convert_to_screen(fract_pos: Vec2) -> Vec2 {
    let (x, y) = ((fract_pos + vec2(1.0, 1.0)) / 2.0).into();
    vec2(x * screen_width(), y * screen_height())
}

/// Convert a velocity in range [-1; 1] to a velocity in pixels.
fn convert_vel_to_screen(fract_vel: Vec2) -> Vec2 {
    let (vx, vy) = (fract_vel / 2.0).into();
    vec2(vx * screen_width(), vy * screen_height())
}

fn handle_exit() -> bool {
    is_key_pressed(KeyCode::Escape) || is_key_pressed(KeyCode::Q)
}

pub fn get_2_mut<T>(v: &mut [T], start: usize, end: usize) -> (&mut T, &mut T) {
    assert_ne!(start, end);
    if end > start {
        let (before, after) = v.split_at_mut(end);
        let c1 = &mut before[start];
        let c2 = &mut after[0];
        (c1, c2)
    } else {
        let (before, after) = v.split_at_mut(start);
        let c1 = &mut before[end];
        let c2 = &mut after[0];
        (c2, c1)
    }
}

pub fn permutations(len: usize) -> impl Iterator<Item = (usize, usize)> {
    (0..len).flat_map(move |i| ((i + 1)..len).map(move |j| (i, j)))
}

fn cfg() -> Conf {
    Conf {
        fullscreen: false,
        window_width: 1280,
        window_height: 720,
        ..Default::default()
    }
}

#[test]
fn test_convert_to_screen() {
    let mock = |fract_pos: Vec2| {
        let (x, y) = ((fract_pos + vec2(1.0, 1.0)) / 2.0).into();
        vec2(x * 1280.0, y * 720.0)
    };
    assert!(mock(vec2(0.0, 0.0)).abs_diff_eq(vec2(640.0, 360.0), 0.0000001));

    let mock = |fract_vel: Vec2| {
        let (vx, vy) = (fract_vel / 2.0).into();
        vec2(vx * 1280.0, vy * 720.0)
    };
    assert!(mock(vec2(0.0, 0.0)).abs_diff_eq(vec2(0.0, 0.0), 0.0000001));
}

#[test]
fn test_permutations() {
    dbg!(permutations(0).collect::<Vec<_>>());
    dbg!(permutations(1).collect::<Vec<_>>());
    dbg!(permutations(2).collect::<Vec<_>>());
    dbg!(permutations(3).collect::<Vec<_>>());
    dbg!(permutations(4).collect::<Vec<_>>());
}
