use rand::distributions::Distribution;
use std::borrow::BorrowMut;
use std::cmp::Ordering;
use std::f64::consts::PI;
use std::sync::{Arc, Mutex};

use crate::raytrace::ray::Ray;
use crate::raytrace::vec3::{dot, Point3, Vec3};

use super::material::Material;

#[derive(Debug, Default, Clone)]
pub struct AABB {
    pub min: Point3,
    pub max: Point3,
}

impl AABB {
    pub fn new(min: Point3, max: Point3) -> AABB {
        AABB { min, max }
    }

    pub fn hit(&self, r: &Ray, mut t_min: f64, mut t_max: f64) -> bool {
        for a in 0..3 {
            let inv_d = 1.0 / r.dir[a];
            let mut t0 = (self.min[a] - r.orig[a]) * inv_d;
            let mut t1 = (self.max[a] - r.orig[a]) * inv_d;
            if inv_d < 0.0 {
                std::mem::swap(&mut t0, &mut t1);
            }

            t_min = if t0 > t_min { t0 } else { t_min };
            t_max = if t1 < t_max { t1 } else { t_max };
            if t_max <= t_min {
                return false;
            }
        }
        true
    }
}

fn surrounding_box(box0: &AABB, box1: &AABB) -> AABB {
    let small = Point3::new(
        if box0.min.x < box1.min.x {
            box0.min.x
        } else {
            box1.min.x
        },
        if box0.min.y < box1.min.y {
            box0.min.y
        } else {
            box1.min.y
        },
        if box0.min.z < box1.min.z {
            box0.min.z
        } else {
            box1.min.z
        },
    );
    let big = Point3::new(
        if box0.max.x > box1.max.x {
            box0.max.x
        } else {
            box1.max.x
        },
        if box0.max.y > box1.max.y {
            box0.max.y
        } else {
            box1.max.y
        },
        if box0.max.z > box1.max.z {
            box0.max.z
        } else {
            box1.max.z
        },
    );
    AABB::new(small, big)
}

pub struct HitRecord {
    pub p: Point3,
    pub normal: Vec3,
    pub material: Arc<dyn Material + Send + Sync>,
    pub t: f64,
    pub u: f64,
    pub v: f64,
    pub front_face: bool,
}

impl HitRecord {
    fn new(
        r: &Ray,
        p: Point3,
        outward_normal: Vec3,
        material: Arc<dyn Material + Send + Sync>,
        t: f64,
        u: f64,
        v: f64,
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
            u,
            v,
            front_face,
        }
    }
}

pub trait Hittable {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB>;
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

    pub fn get_sphere_uv(p: Point3) -> (f64, f64) {
        let theta = f64::acos(-p.y);
        let phi = f64::atan2(-p.z, p.x) + PI;

        (phi / (2.0 * PI), theta / PI)
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
        let outward_normal = (r.at(root) - self.center) / self.radius;
        let (u, v) = Sphere::get_sphere_uv(outward_normal);

        Some(HitRecord::new(
            r,
            r.at(root),
            outward_normal,
            self.material.clone(),
            root,
            u,
            v,
        ))
    }

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        Some(AABB::new(
            self.center - Vec3::new(self.radius, self.radius, self.radius),
            self.center + Vec3::new(self.radius, self.radius, self.radius),
        ))
    }
}

pub struct MovingSphere {
    center0: Point3,
    center1: Point3,
    t0: f64,
    t1: f64,
    radius: f64,
    material: Arc<dyn Material + Send + Sync>,
}

impl MovingSphere {
    pub fn new(
        center0: Point3,
        center1: Point3,
        t0: f64,
        t1: f64,
        radius: f64,
        material: Arc<dyn Material + Send + Sync>,
    ) -> MovingSphere {
        MovingSphere {
            center0,
            center1,
            t0,
            t1,
            radius,
            material,
        }
    }

    pub fn center(&self, time: f64) -> Point3 {
        self.center0 + ((time - self.t0) / (self.t1 - self.t0)) * (self.center1 - self.center0)
    }
}

impl Hittable for MovingSphere {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc = r.orig - self.center(r.time);
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
            (r.at(root) - self.center(r.time)) / self.radius,
            self.material.clone(),
            root,
            0.0,
            0.0,
        ))
    }

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        let box0 = AABB::new(
            self.center(t0) - Vec3::new(self.radius, self.radius, self.radius),
            self.center(t0) + Vec3::new(self.radius, self.radius, self.radius),
        );
        let box1 = AABB::new(
            self.center(t1) - Vec3::new(self.radius, self.radius, self.radius),
            self.center(t1) + Vec3::new(self.radius, self.radius, self.radius),
        );

        Some(surrounding_box(&box0, &box1))
    }
}

pub struct BVHNode {
    left: Arc<dyn Hittable + Send + Sync>,
    right: Arc<dyn Hittable + Send + Sync>,
    box_: AABB,
}

impl BVHNode {
    pub fn new(rng: &mut dyn rand::RngCore, world: HittableList, t0: f64, t1: f64) -> BVHNode {
        Self::new_impl(rng, world.objects.to_vec(), t0, t1)
    }

    fn new_impl(
        rng: &mut dyn rand::RngCore,
        mut list: Vec<Arc<dyn Hittable + Send + Sync>>,
        t0: f64,
        t1: f64,
    ) -> BVHNode {
        let uni = rand::distributions::Uniform::from(0..3);
        let axis = uni.sample(rng);
        let comparator = |a: &Arc<dyn Hittable + Send + Sync>,
                          b: &Arc<dyn Hittable + Send + Sync>| {
            let box_a = a.bounding_box(t0, t1).unwrap();
            let box_b = b.bounding_box(t0, t1).unwrap();
            match axis {
                0 => box_a.min.x.partial_cmp(&box_b.min.x).unwrap(),
                1 => box_a.min.y.partial_cmp(&box_b.min.y).unwrap(),
                2 => box_a.min.z.partial_cmp(&box_b.min.z).unwrap(),
                _ => panic!("invalid axis"),
            }
        };

        let (left, right) = if list.len() == 1 {
            (list[0].clone(), list[0].clone())
        } else if list.len() == 2 {
            if comparator(&list[0], &list[1]) == Ordering::Less {
                (list[0].clone(), list[1].clone())
            } else {
                (list[1].clone(), list[0].clone())
            }
        } else {
            list.sort_by(comparator);
            let mid = list.len() / 2;
            (
                Arc::new(BVHNode::new_impl(rng, list[..mid].to_vec(), t0, t1))
                    as Arc<dyn Hittable + Send + Sync>,
                Arc::new(BVHNode::new_impl(rng, list[mid..].to_vec(), t0, t1))
                    as Arc<dyn Hittable + Send + Sync>,
            )
        };

        let box_left = left.bounding_box(t0, t1).unwrap();
        let box_right = right.bounding_box(t0, t1).unwrap();
        let box_ = surrounding_box(&box_left, &box_right);

        BVHNode { left, right, box_ }
    }
}

impl Hittable for BVHNode {
    fn hit(&self, r: &Ray, t_min: f64, mut t_max: f64) -> Option<HitRecord> {
        if !self.box_.hit(r, t_min, t_max) {
            None
        } else {
            let left_rec = self.left.hit(r, t_min, t_max);
            if let Some(rec) = &left_rec {
                t_max = f64::min(t_max, rec.t);
            }

            self.right.hit(r, t_min, t_max).or_else(|| left_rec)
        }
    }
    fn bounding_box(&self, _t0: f64, _t1: f64) -> Option<AABB> {
        Some(self.box_.clone())
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

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        self.objects
            .iter()
            .filter_map(|object| object.bounding_box(t0, t1))
            .reduce(|a, b| surrounding_box(&a, &b))
    }
}
