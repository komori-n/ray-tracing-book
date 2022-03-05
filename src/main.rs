const WIDTH: usize = 256;
const HEIGHT: usize = 256;

fn output() {
    println!("P3");
    println!("{} {}", WIDTH, HEIGHT);
    println!("255");

    for j in (0..HEIGHT).rev() {
        eprint!("\rScanlines remaining: {}   ", j);
        for i in 0..WIDTH {
            let r = (i as f64) / ((WIDTH - 1) as f64);
            let g = (j as f64) / ((HEIGHT - 1) as f64);
            let b = 0.25;

            let ir = (r * 255.999) as i64;
            let ig = (g * 255.999) as i64;
            let ib = (b * 255.999) as i64;
            println!("{} {} {}", ir, ig, ib);
        }
    }
    eprintln!("");
    eprintln!("Done");
}

fn main() {
    output()
}
