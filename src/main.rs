mod raytrace;

use std::intrinsics::discriminant_value;

use raytrace::vec3::unit;

use crate::raytrace::color::{ppm_string, Color};
use crate::raytrace::ray::Ray;
use crate::raytrace::vec3::{dot, Point3, Vec3};

const ASPECT_RATIO: f64 = 16.0 / 9.0;
const WIDTH: usize = 400;
const HEIGHT: usize = ((WIDTH as f64) / ASPECT_RATIO) as usize;

fn output() {
    let viewpoint_height = 2.0;
    let viewpoint_width = ASPECT_RATIO * viewpoint_height;
    let focal_length = 1.0;

    let origin = Point3::new(0.0, 0.0, 0.0);
    let horizontal = Vec3::new(viewpoint_width, 0.0, 0.0);
    let vertical = Vec3::new(0.0, viewpoint_height, 0.0);
    let lower_left_corner =
        origin - horizontal / 2.0 - vertical / 2.0 - Vec3::new(0.0, 0.0, focal_length);

    println!("P3");
    println!("{} {}", WIDTH, HEIGHT);
    println!("255");

    for j in (0..HEIGHT).rev() {
        eprint!("\rScanlines remaining: {}   ", j);
        for i in 0..WIDTH {
            let u = (i as f64) / ((WIDTH - 1) as f64);
            let v = (j as f64) / ((HEIGHT - 1) as f64);
            let ray = Ray::new(
                origin,
                lower_left_corner + u * horizontal + v * vertical - origin,
            );
            let c = ray_color(&ray);
            println!("{}", ppm_string(c));
        }
    }
    eprintln!("");
    eprintln!("Done");
}

fn ray_color(ray: &Ray) -> Color {
    let t = hit_sphere(&Point3::new(0.0, 0.0, -1.0), 0.5, ray);
    if t > 0.0 {
        let n = unit(ray.at(t) - Vec3::new(0.0, 0.0, -1.0));
        return 0.5 * Color::new(n.x + 1.0, n.y + 1.0, n.z + 1.0);
    }
    let unit_direction = unit(ray.dir);
    let t = 0.5 * (unit_direction.y + 1.0);

    (1.0 - t) * Color::new(1.0, 1.0, 1.0) + t * Color::new(0.5, 0.7, 1.0)
}

fn hit_sphere(center: &Point3, radius: f64, r: &Ray) -> f64 {
    let oc = r.orig - *center;
    let a = dot(r.dir, r.dir);
    let half_b = dot(oc, r.dir);
    let c = dot(oc, oc) - radius * radius;
    let discriminant = half_b * half_b - a * c;

    if discriminant < 0.0 {
        -1.0
    } else {
        (-half_b - f64::sqrt(discriminant)) / a
    }
}

fn main() {
    output();
}
