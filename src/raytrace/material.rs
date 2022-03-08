use std::sync::Arc;

use crate::raytrace::{
    color::Color,
    hittable::HitRecord,
    ray::Ray,
    vec3::{dot, random_in_unit_sphere, random_unit_vector, reflect, refract, unit},
};

use super::{
    texture::{SolidColor, Texture},
    vec3::Point3,
};

pub trait Material: Sync + Send {
    fn scatter(
        &self,
        rng: &mut dyn rand::RngCore,
        r_in: &Ray,
        rec: &HitRecord,
    ) -> Option<(Ray, Color)>;

    fn emitted(&self, u: f64, v: f64, p: &Point3) -> Color {
        Color::default()
    }
}

pub struct Lambertian {
    albedo: Arc<dyn Texture + Send + Sync>,
}

impl Lambertian {
    pub fn new(albedo: Color) -> Lambertian {
        Lambertian {
            albedo: Arc::new(SolidColor::new(albedo)),
        }
    }

    pub fn from_texture(a: Arc<dyn Texture + Send + Sync>) -> Lambertian {
        Lambertian { albedo: a }
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

        let scattered = Ray::new(rec.p, scatter_direction, r_in.time);
        let attenuation = self.albedo.value(rec.u, rec.v, &rec.p);
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
        let scattered = Ray::new(
            rec.p,
            reflected + self.fuzz * random_in_unit_sphere(rng),
            r_in.time,
        );
        let attenuation = self.albedo;

        if dot(scattered.dir, rec.normal) > 0.0 {
            Some((scattered, attenuation))
        } else {
            None
        }
    }
}

pub struct Dielectric {
    ir: f64,
}

impl Dielectric {
    pub fn new(index_of_refraction: f64) -> Dielectric {
        Dielectric {
            ir: index_of_refraction,
        }
    }
}

impl Material for Dielectric {
    fn scatter(
        &self,
        rng: &mut dyn rand::RngCore,
        r_in: &Ray,
        rec: &HitRecord,
    ) -> Option<(Ray, Color)> {
        let attenuation = Color::new(1.0, 1.0, 1.0);
        let refraction_ratio = if rec.front_face {
            1.0 / self.ir
        } else {
            self.ir
        };

        let unit_direction = unit(r_in.dir);
        let cos_theta = f64::min(dot(-unit_direction, rec.normal), 1.0);
        let sin_theta = f64::sqrt(1.0 - cos_theta * cos_theta);

        let direction = if refraction_ratio * sin_theta > 1.0 {
            reflect(unit_direction, rec.normal)
        } else {
            refract(unit_direction, rec.normal, refraction_ratio)
        };

        let scattered = Ray::new(rec.p, direction, r_in.time);
        Some((scattered, attenuation))
    }
}

pub struct DiffuseLight {
    emit: Arc<dyn Texture + Send + Sync>,
}

impl DiffuseLight {
    pub fn new(emit: Color) -> DiffuseLight {
        DiffuseLight {
            emit: Arc::new(SolidColor::new(emit)),
        }
    }

    pub fn from_texture(a: Arc<dyn Texture + Send + Sync>) -> DiffuseLight {
        DiffuseLight { emit: a }
    }
}

impl Material for DiffuseLight {
    fn scatter(
        &self,
        _rng: &mut dyn rand::RngCore,
        _r_in: &Ray,
        _rec: &HitRecord,
    ) -> Option<(Ray, Color)> {
        None
    }

    fn emitted(&self, u: f64, v: f64, p: &Point3) -> Color {
        self.emit.value(u, v, p)
    }
}
