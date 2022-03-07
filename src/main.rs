mod raytrace;

use std::f64::consts::PI;
use std::mem::Discriminant;
use std::rc::Rc;
use std::sync::{Arc, Mutex};

use itertools::iproduct;
use rand::prelude::*;
use rayon::prelude::*;

use raytrace::hittable::{Hittable, MovingSphere};
use raytrace::vec3::unit;

use crate::raytrace::camera::Camera;
use crate::raytrace::color::{ppm_string, Color};
use crate::raytrace::hittable::{HittableList, Sphere};
use crate::raytrace::material::{Dielectric, Lambertian, Material, Metal};
use crate::raytrace::ray::Ray;
use crate::raytrace::vec3::{Point3, Vec3};

const ASPECT_RATIO: f64 = 16.0 / 9.0;
const WIDTH: usize = 400;
const HEIGHT: usize = ((WIDTH as f64) / ASPECT_RATIO) as usize;
const SAMPLES_PER_PIXEL: i64 = 100;
const MAX_DEPTH: i64 = 50;
const APERTURE: f64 = 0.1;

fn random_scene(rng: &mut impl rand::RngCore) -> HittableList {
    let mut objects: Vec<Arc<dyn Hittable + Send + Sync>> = Vec::new();
    let ground_material: Arc<dyn Material + Send + Sync> =
        Arc::new(Lambertian::new(Color::new(0.5, 0.5, 0.5)));
    let sphere = Sphere::new(Point3::new(0.0, -1000.0, 0.0), 1000.0, ground_material);
    objects.push(Arc::new(sphere));

    let uni = rand::distributions::Uniform::from(0.0..1.0);
    for a in -11..11 {
        for b in -11..11 {
            let choose_mat = uni.sample(rng);
            let center = Point3::new(
                (a as f64) + 0.9 * uni.sample(rng),
                0.2,
                (b as f64) + 0.9 * uni.sample(rng),
            );

            if (center - Point3::new(4.0, 0.2, 0.0)).length() > 0.9 {
                let material: Arc<dyn Material + Send + Sync> = if choose_mat < 0.8 {
                    let albedo = Color::random(rng) * Color::random(rng);
                    Arc::new(Lambertian::new(albedo))
                } else if choose_mat < 0.95 {
                    let albedo = Color::random_with_minmax(rng, 0.5, 1.0);
                    let fuzz = uni.sample(rng) / 0.5;
                    Arc::new(Metal::new(albedo, fuzz))
                } else {
                    Arc::new(Dielectric::new(1.5))
                };

                if choose_mat < 0.8 {
                    let center2 = center + Vec3::new(0.0, 0.5 * uni.sample(rng), 0.0);
                    objects.push(Arc::new(MovingSphere::new(
                        center, center2, 0.0, 1.0, 0.2, material,
                    )));
                } else {
                    objects.push(Arc::new(Sphere::new(center, 0.2, material)));
                }
            }
        }
    }

    let material1: Arc<dyn Material + Send + Sync> = Arc::new(Dielectric::new(1.5));
    objects.push(Arc::new(Sphere::new(
        Point3::new(0.0, 1.0, 0.0),
        1.0,
        material1,
    )));

    let material2: Arc<dyn Material + Send + Sync> =
        Arc::new(Lambertian::new(Color::new(0.4, 0.2, 0.1)));
    objects.push(Arc::new(Sphere::new(
        Point3::new(-4.0, 1.0, 0.0),
        1.0,
        material2,
    )));

    let material3: Arc<dyn Material + Send + Sync> =
        Arc::new(Metal::new(Color::new(0.7, 0.6, 0.5), 0.0));
    objects.push(Arc::new(Sphere::new(
        Point3::new(4.0, 1.0, 0.0),
        1.0,
        material3,
    )));

    HittableList::new(objects)
}

fn output() {
    let mut rng = rand::thread_rng();
    let uni = rand::distributions::Uniform::from(0.0..1.0);

    let mut world = random_scene(&mut rng);

    let lookfrom = Point3::new(13.0, 2.0, 3.0);
    let lookat = Point3::new(0.0, 0.0, 0.0);
    let vup = Vec3::new(0.0, 1.0, 0.0);
    let dist_to_focus = 10.0;
    let cam = Camera::new(
        lookfrom,
        lookat,
        vup,
        20.0,
        ASPECT_RATIO,
        APERTURE,
        dist_to_focus,
        0.0,
        1.0,
    );

    println!("P3");
    println!("{} {}", WIDTH, HEIGHT);
    println!("255");

    let a: Vec<_> = iproduct!((0..HEIGHT).rev(), 0..WIDTH).collect();
    let x: Vec<_> = a
        .into_par_iter()
        .map(|(j, i)| {
            let mut rng = rand::thread_rng();
            let mut color = Color::default();
            for _ in 0..SAMPLES_PER_PIXEL {
                let u = ((i as f64) + uni.sample(&mut rng)) / ((WIDTH - 1) as f64);
                let v = ((j as f64) + uni.sample(&mut rng)) / ((HEIGHT - 1) as f64);
                let r = cam.get_ray(&mut rng, u, v);
                color += ray_color(&mut rng, &r, &world, MAX_DEPTH);
            }
            color
        })
        .collect();
    x.iter()
        .for_each(|c| println!("{}", ppm_string(*c, SAMPLES_PER_PIXEL)));
    eprintln!("");
    eprintln!("Done");
}

fn ray_color(rng: &mut dyn rand::RngCore, ray: &Ray, world: &dyn Hittable, depth: i64) -> Color {
    if depth <= 0 {
        return Color::default();
    }

    if let Some(rec) = world.hit(ray, 0.001, f64::MAX) {
        if let Some((scattered, attenuation)) = rec.material.scatter(rng, ray, &rec) {
            return attenuation * ray_color(rng, &scattered, world, depth - 1);
        }
        return Color::default();
    }
    let unit_direction = unit(ray.dir);
    let t = 0.5 * (unit_direction.y + 1.0);
    return (1.0 - t) * Color::new(1.0, 1.0, 1.0) + t * Color::new(0.5, 0.7, 1.0);
}

fn main() {
    output();
}
