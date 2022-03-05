use crate::raytrace::vec3::Vec3;

pub type Color = Vec3;

pub fn ppm_string(c: Color) -> String {
    let ir = (c.x * 255.999) as i64;
    let ig = (c.y * 255.999) as i64;
    let ib = (c.z * 255.999) as i64;

    format!("{} {} {}", ir, ig, ib)
}
