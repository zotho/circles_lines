use macroquad::prelude::*;

pub fn get_center(points: &[&mut Vec2]) -> Vec2 {
    points.iter().map(|p| &**p).sum::<Vec2>() / points.len() as f32
}

pub fn set_center(points: &mut [&mut Vec2], center: Vec2) {
    let actual = get_center(&points);
    let diff = center - actual; // A -> C
    points.iter_mut().for_each(|p| **p += diff);
}

#[test]
fn test_get_center() {
    assert!(Vec2::abs_diff_eq(
        get_center(&vec![&mut vec2(0.0, 1.0), &mut vec2(1.0, 0.0)]),
        vec2(0.5, 0.5),
        std::f32::EPSILON
    ));
}

#[test]
fn test_set_center() {
    let mut a = vec2(0.0, 1.0);
    let mut b = vec2(1.0, 0.0);
    let mut points = vec![&mut a, &mut b];
    set_center(&mut points, vec2(0.0, 0.0));
    assert!(Vec2::abs_diff_eq(a, vec2(-0.5, 0.5), std::f32::EPSILON));
    assert!(Vec2::abs_diff_eq(b, vec2(0.5, -0.5), std::f32::EPSILON));
}
