mod raytrace;

use rand::prelude::*;

use raytrace::hittable::Hittable;
use raytrace::vec3::{random_in_hemisphere, random_in_unit_sphere, random_unit_veector, unit};

use crate::raytrace::camera::Camera;
use crate::raytrace::color::{ppm_string, Color};
use crate::raytrace::hittable::{HittableList, Sphere};
use crate::raytrace::ray::Ray;
use crate::raytrace::vec3::Point3;

const ASPECT_RATIO: f64 = 16.0 / 9.0;
const WIDTH: usize = 400;
const HEIGHT: usize = ((WIDTH as f64) / ASPECT_RATIO) as usize;
const SAMPLES_PER_PIXEL: i64 = 100;
const MAX_DEPTH: i64 = 50;

fn output() {
    let mut rng = rand::thread_rng();
    let uni = rand::distributions::Uniform::from(0.0..1.0);

    let mut world = HittableList::default();
    world.add(Box::new(Sphere::new(Point3::new(0.0, 0.0, -1.0), 0.5)));
    world.add(Box::new(Sphere::new(Point3::new(0.0, -100.5, -1.0), 100.0)));

    let cam = Camera::new();

    println!("P3");
    println!("{} {}", WIDTH, HEIGHT);
    println!("255");

    for j in (0..HEIGHT).rev() {
        eprint!("\rScanlines remaining: {}   ", j);
        for i in 0..WIDTH {
            let mut color = Color::default();
            for _ in 0..SAMPLES_PER_PIXEL {
                let u = ((i as f64) + uni.sample(&mut rng)) / ((WIDTH - 1) as f64);
                let v = ((j as f64) + uni.sample(&mut rng)) / ((HEIGHT - 1) as f64);
                let r = cam.get_ray(u, v);
                color += ray_color(&mut rng, &r, &world, MAX_DEPTH);
            }
            println!("{}", ppm_string(color, SAMPLES_PER_PIXEL));
        }
    }
    eprintln!("");
    eprintln!("Done");
}

fn ray_color(rng: &mut dyn rand::RngCore, ray: &Ray, world: &dyn Hittable, depth: i64) -> Color {
    if depth <= 0 {
        return Color::default();
    }

    if let Some(rec) = world.hit(ray, 0.001, f64::MAX) {
        let target = rec.p + rec.normal + random_in_hemisphere(rng, &rec.normal);
        return 0.5 * ray_color(rng, &Ray::new(rec.p, target - rec.p), world, depth - 1);
    }
    let unit_direction = unit(ray.dir);
    let t = 0.5 * (unit_direction.y + 1.0);
    return (1.0 - t) * Color::new(1.0, 1.0, 1.0) + t * Color::new(0.5, 0.7, 1.0);
}

fn main() {
    output();
}
