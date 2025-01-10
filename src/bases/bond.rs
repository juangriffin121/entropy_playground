use super::{atom::Atom, physics::Force, recipient::Recipient};
use serde::Deserialize;
use std::{collections::HashMap, fs};

#[derive(Debug)]
pub struct Bond {
    pub atom1_id: usize,
    pub atom2_id: usize,
    pub k: f32, // Spring constant
    pub equilibrium_distance: f32,
    pub force: (f32, f32), //cached force value
    pub breaking_distance: f32,
}

pub struct BondForce {}

impl Bond {
    pub fn new(
        atom1_id: usize,
        atom2_id: usize,
        k: f32,
        equilibrium_distance: f32,
        breaking_distance: f32,
    ) -> Self {
        Bond {
            atom1_id,
            atom2_id,
            k,
            equilibrium_distance,
            force: (0.0, 0.0),
            breaking_distance,
        }
    }

    pub fn is_atom1(&self, atom: &Atom) -> bool {
        if atom.id == self.atom1_id {
            return true;
        } else if atom.id == self.atom2_id {
            return false;
        } else {
            unreachable!()
        }
    }
}

impl Force for BondForce {
    fn calculate(&self, atom: &Atom, recipient: &Recipient) -> (f32, f32) {
        let mut fx = 0.0;
        let mut fy = 0.0;
        for bond_idx in &atom.bonds {
            let bond = &recipient.bond_registry[&bond_idx];
            if bond.is_atom1(atom) {
                let f = bond.force;
                fx += f.0;
                fy += f.1;
            } else {
                let f = bond.force;
                fx -= f.0;
                fy -= f.1;
            }
        }
        (fx, fy)
    }
}

#[derive(Debug, Deserialize, Clone)]
pub struct BondProperties {
    pub strength: f32, // Strength of the bond
    pub breaking_distance: f32,
}

#[derive(Debug, Deserialize)]
pub struct BondDefinitions {
    pub data: HashMap<String, BondProperties>,
}

impl BondDefinitions {
    pub fn load(path: &str) -> Self {
        let json_content = fs::read_to_string(path).expect("Failed to read JSON file");
        serde_json::from_str(&json_content).expect("Failed to parse JSON")
    }

    pub fn to_tuple_map(&self) -> HashMap<(String, String), BondProperties> {
        let mut tuple_map = HashMap::new();

        for (key, properties) in &self.data {
            let parts: Vec<&str> = key.split('-').collect();

            if parts.len() == 2 {
                let (element1, element2) = if parts[0] <= parts[1] {
                    (parts[0], parts[1])
                } else {
                    (parts[1], parts[0])
                };

                tuple_map.insert(
                    (element1.to_string(), element2.to_string()),
                    properties.clone(),
                );
                tuple_map.insert(
                    (parts[1].to_string(), parts[0].to_string()),
                    properties.clone(),
                );
            } else {
                eprintln!("Invalid bond format: {}", key);
            }
        }
        tuple_map
    }
}
