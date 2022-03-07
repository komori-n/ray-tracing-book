use std::borrow::BorrowMut;
use std::sync::{Arc, Mutex};

use crate::raytrace::ray::Ray;
use crate::raytrace::vec3::{dot, Point3, Vec3};

use super::material::Material;

pub struct HitRecord {
    pub p: Point3,
    pub normal: Vec3,
    pub material: Arc<dyn Material + Send + Sync>,
    pub t: f64,
    pub front_face: bool,
}

impl HitRecord {
    fn new(
        r: &Ray,
        p: Point3,
        outward_normal: Vec3,
        material: Arc<dyn Material + Send + Sync>,
        t: f64,
    ) -> HitRecord {
        let front_face = dot(r.dir, outward_normal) < 0.0;
        let normal = if front_face {
            outward_normal
        } else {
            -outward_normal
        };

        HitRecord {
            p,
            normal,
            material,
            t,
            front_face,
        }
    }
}

pub trait Hittable {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
}

pub struct Sphere {
    center: Point3,
    radius: f64,
    material: Arc<dyn Material + Send + Sync>,
}

impl Sphere {
    pub fn new(center: Point3, radius: f64, material: Arc<dyn Material + Send + Sync>) -> Sphere {
        Sphere {
            center,
            radius,
            material,
        }
    }
}

impl Hittable for Sphere {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc = r.orig - self.center;
        let a = r.dir.length_squared();
        let half_b = dot(oc, r.dir);
        let c = oc.length_squared() - self.radius * self.radius;

        let discriminant = half_b * half_b - a * c;
        if discriminant < 0.0 {
            return None;
        }

        let sqrtd = f64::sqrt(discriminant);
        let mut root = (-half_b - sqrtd) / a;
        if root < t_min || t_max < root {
            root = (-half_b + sqrtd) / a;
            if root < t_min || t_max < root {
                return None;
            }
        }
        let root = root;

        Some(HitRecord::new(
            r,
            r.at(root),
            (r.at(root) - self.center) / self.radius,
            self.material.clone(),
            root,
        ))
    }
}

#[derive(Default)]
pub struct HittableList {
    objects: Arc<Vec<Arc<dyn Hittable + Send + Sync>>>,
}

impl HittableList {
    pub fn new(objects: Vec<Arc<dyn Hittable + Send + Sync>>) -> HittableList {
        HittableList {
            objects: Arc::new(objects),
        }
    }
}

impl Hittable for HittableList {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut closest_so_far = t_max;
        let mut hit_record = None;

        for object in self.objects.iter() {
            if let Some(hit) = object.hit(r, t_min, closest_so_far) {
                closest_so_far = hit.t;
                hit_record = Some(hit);
            }
        }

        hit_record
    }
}
