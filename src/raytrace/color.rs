use crate::raytrace::vec3::Vec3;

pub type Color = Vec3;

fn clamp(input: f64, min: f64, max: f64) -> f64 {
    if input < min {
        min
    } else if input > max {
        max
    } else {
        input
    }
}

pub fn ppm_string(c: Color, samples_per_pixel: i64) -> String {
    let r = c.x / (samples_per_pixel as f64);
    let g = c.y / (samples_per_pixel as f64);
    let b = c.z / (samples_per_pixel as f64);

    let ir = (256.0 * clamp(r, 0.0, 0.999)) as i64;
    let ig = (256.0 * clamp(g, 0.0, 0.999)) as i64;
    let ib = (256.0 * clamp(b, 0.0, 0.999)) as i64;

    format!("{} {} {}", ir, ig, ib)
}
