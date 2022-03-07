use crate::raytrace::{
    color::Color,
    hittable::HitRecord,
    ray::Ray,
    vec3::{dot, random_unit_vector, reflect, unit},
};

use super::vec3::random_in_unit_sphere;

pub trait Material {
    fn scatter(
        &self,
        rng: &mut dyn rand::RngCore,
        r_in: &Ray,
        rec: &HitRecord,
    ) -> Option<(Ray, Color)>;
}

pub struct Lambertian {
    albedo: Color,
}

impl Lambertian {
    pub fn new(albedo: Color) -> Lambertian {
        Lambertian { albedo }
    }
}

impl Material for Lambertian {
    fn scatter(
        &self,
        rng: &mut dyn rand::RngCore,
        r_in: &Ray,
        rec: &HitRecord,
    ) -> Option<(Ray, Color)> {
        let mut scatter_direction = rec.normal + random_unit_vector(rng);
        if scatter_direction.near_zero() {
            scatter_direction = rec.normal;
        }

        let scattered = Ray::new(rec.p, scatter_direction);
        let attenuation = self.albedo;
        Some((scattered, attenuation))
    }
}

pub struct Metal {
    albedo: Color,
    fuzz: f64,
}

impl Metal {
    pub fn new(albedo: Color, fuzz: f64) -> Metal {
        Metal {
            albedo,
            fuzz: if fuzz < 1.0 { fuzz } else { 1.0 },
        }
    }
}

impl Material for Metal {
    fn scatter(
        &self,
        rng: &mut dyn rand::RngCore,
        r_in: &Ray,
        rec: &HitRecord,
    ) -> Option<(Ray, Color)> {
        let reflected = reflect(unit(r_in.dir), rec.normal);
        let scattered = Ray::new(rec.p, reflected + self.fuzz * random_in_unit_sphere(rng));
        let attenuation = self.albedo;

        if dot(scattered.dir, rec.normal) > 0.0 {
            Some((scattered, attenuation))
        } else {
            None
        }
    }
}
