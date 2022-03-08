use std::sync::Arc;

use crate::raytrace::{color::Color, vec3::Point3};

use super::perlin::Perlin;

pub trait Texture {
    fn value(&self, u: f64, v: f64, p: &Point3) -> Color;
}

#[derive(Default)]
pub struct SolidColor {
    pub color: Color,
}

impl SolidColor {
    pub fn new(color: Color) -> SolidColor {
        SolidColor { color }
    }
}

impl Texture for SolidColor {
    fn value(&self, u: f64, v: f64, p: &Point3) -> Color {
        self.color
    }
}

pub struct CheckerTexture {
    pub odd: Arc<dyn Texture + Send + Sync>,
    pub even: Arc<dyn Texture + Send + Sync>,
}

impl CheckerTexture {
    pub fn new(even: Color, odd: Color) -> CheckerTexture {
        CheckerTexture {
            odd: Arc::new(SolidColor::new(odd)),
            even: Arc::new(SolidColor::new(even)),
        }
    }

    pub fn from_texture(
        even: Arc<dyn Texture + Send + Sync>,
        odd: Arc<dyn Texture + Send + Sync>,
    ) -> CheckerTexture {
        CheckerTexture { odd, even }
    }
}

impl Texture for CheckerTexture {
    fn value(&self, u: f64, v: f64, p: &Point3) -> Color {
        let sines = (10.0 * p.x).sin() * (10.0 * p.y).sin() * (10.0 * p.z).sin();
        if sines < 0.0 {
            self.odd.value(u, v, p)
        } else {
            self.even.value(u, v, p)
        }
    }
}

pub struct NoiseTexture {
    noise: Perlin,
    scale: f64,
}

impl NoiseTexture {
    pub fn new(rng: &mut dyn rand::RngCore, scale: f64) -> NoiseTexture {
        NoiseTexture {
            noise: Perlin::new(rng),
            scale,
        }
    }
}

impl Texture for NoiseTexture {
    fn value(&self, u: f64, v: f64, p: &Point3) -> Color {
        Color::new(1.0, 1.0, 1.0)
            * 0.5
            * (1.0 + (self.scale * p.z + 10.0 * self.noise.turb(p, 7)).sin())
    }
}
