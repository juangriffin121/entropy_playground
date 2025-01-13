use super::molecule::Molecule;
use super::recipient::Recipient;
use plotters::style::RGBColor;
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::fs;

#[derive(Debug)]
pub struct Atom {
    pub id: usize,
    pub mass: f32,
    pub radius: f32,
    pub color: RGBColor,
    pub position: (f32, f32),
    pub velocity: (f32, f32),
    pub molecule_id: Option<usize>,
    pub id_in_molecule: Option<usize>,
    pub chemical_nature: String,
    pub bonded_atoms: HashSet<usize>,
}

#[derive(Debug, Deserialize, Clone)]
pub struct Element {
    pub color: [u8; 3],
    pub mass: f32,
    pub radius: f32,
}

#[derive(Debug, Deserialize)]
pub struct ElementsJSON {
    pub data: HashMap<String, Element>,
}

impl ElementsJSON {
    pub fn load(path: &str) -> Self {
        let json_content = fs::read_to_string(path).expect("Failed to read JSON file");
        serde_json::from_str(&json_content).expect("Failed to parse JSON")
    }
}

impl Atom {
    pub fn new(
        id: usize,
        position: (f32, f32),
        velocity: (f32, f32),
        element: &str,
        elements_data: &HashMap<String, Element>,
    ) -> Self {
        let element_data = &elements_data[element];
        let (mass, radius, [r, g, b]) =
            (element_data.mass, element_data.radius, element_data.color);
        let color = RGBColor(r, g, b);
        let molecule_id = None;
        let id_in_molecule = None;
        let chemical_nature = element.to_string();
        let bonded_atoms = HashSet::new();

        Self {
            id,
            position,
            velocity,
            mass,
            radius,
            color,
            molecule_id,
            id_in_molecule,
            chemical_nature,
            bonded_atoms,
        }
    }
}

#[derive(Debug, Deserialize)]
pub struct DistributionJSON {
    pub map: HashMap<String, i32>,
}

impl DistributionJSON {
    pub fn load(path: &str) -> Self {
        let json_content = fs::read_to_string(path).expect("Failed to read JSON file");
        serde_json::from_str(&json_content).expect("Failed to parse JSON")
    }
    pub fn to_cdf(&self) -> Vec<(String, f32)> {
        let mut cdf = Vec::new();
        let mut sum = 0.0;
        let total: f32 = self.map.values().map(|&weight| weight as f32).sum();
        for (name, prob) in self.map.iter() {
            sum += *prob as f32 / total;
            cdf.push((name.clone(), sum));
        }
        cdf.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        cdf
    }
}
