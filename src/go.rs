use std::fmt::Debug;
use std::fmt::Formatter;
use std::fmt::Error;
use std::fmt::Display;
use std::collections::HashSet;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Stone {
    Black,
    White,
}

impl Stone {
    fn is_friend(&self, other: &Option<Stone>) -> bool {
        Some(*self) == *other
    }

    fn is_enemy(&self, other: &Option<Stone>) -> bool {
        match other {
            Some(o) => self != o,
            None => false
        }
    }
}

impl Debug for Stone {
    fn fmt(&self, f: &mut Formatter) -> Result<(), Error> {
        write!(f, "{}", self)
    }
}

impl Display for Stone {
    fn fmt(&self, f: &mut Formatter) -> Result<(), Error> {
        write!(f, "{}", match self {
            Stone::Black => 'B',
            Stone::White => 'W',
        })
    }
}

struct StandardBoard([Option<Stone>; 361]);

impl Debug for StandardBoard {
    fn fmt(&self, f: &mut Formatter) -> Result<(), Error> {
        for r in 0..19 {
            for c in 0..19 {
                match self.0[self.index(r, c)] {
                    Some(s) => write!(f, "{}", s)?,
                    None => write!(f, "+")?,
                }
            }
            writeln!(f, "")?;
        }
        Ok(())
    }
}

impl StandardBoard {
    fn new() -> StandardBoard {
        StandardBoard([None; 361])
    }

    fn index(&self, row: usize, col: usize) -> usize { row * 19 + col }

    fn get_stone(&self, row: usize, col: usize) -> Option<Stone> {
        self.0[self.index(row, col)]
    }

    fn get(&self, row: usize, col: usize) -> Intersection {
        Intersection::new(&self, row, col)
    }

    fn set(&mut self, row: usize, col: usize, stone: Option<Stone>) {
        self.0[self.index(row, col)] = stone;
    }
}

#[derive(Debug, Clone)]
struct Intersection<'a> {
    board: &'a StandardBoard,
    row: usize,
    col: usize,
}

impl <'a>Intersection<'a> {
    fn new(board: &'a StandardBoard, row: usize, col: usize) -> Intersection<'a> {
        Intersection { board, row, col }
    }

    fn stone(&self) -> Option<Stone> {
        self.board.get_stone(self.row, self.col)
    }

    fn neighbors(&self) -> Vec<Intersection> {
        let indices = match (self.row, self.col) {
            (0, 0) => vec![(0, 1), (1, 0)],
            (0, 18) => vec![(0, 17), (1, 18)],
            (18, 0) => vec![(17, 0), (18, 1)],
            (18, 18) => vec![(17, 18), (18, 17)],
            (0, c) => vec![(0, c - 1), (0, c + 1), (1, c)],
            (18, c) => vec![(18, c - 1), (18, c + 1), (17, c)],
            (r, 0) => vec![(r - 1, 0), (r + 1, 0), (r, 1)],
            (r, 18) => vec![(r - 1, 18), (r + 1, 18), (r, 17)],
            (r, c) => vec![(r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)]
        };
        indices.into_iter().map(|(r, c)| self.board.get(r, c)).collect()
    }

    fn is_alive(&'a self) -> bool {
        match self.stone() {
            None => true, // empty intersections are (vacuously) alive
            Some(s) => {
                let neighbors = self.neighbors();
                // any empty neighbors means the stone is alive
                if neighbors.iter().any(|n| n.stone().is_none()) {
                    return true
                }
                // all enemy neighbors means the stone is dead
                if neighbors.iter().all(|n| {s.is_enemy(&n.stone())}) {
                    return false
                }
                // if no empty neighbors, but not all enemies, check friendly neighbors for liveness
                let mut visited = HashSet::new();
                let mut to_visit = vec![(self.row, self.col)];
                while !to_visit.is_empty() {
                    let (r, c) = to_visit.pop().unwrap();
                    let i = self.board.get(r, c);
                    if visited.insert((i.row, i.col)) {
                        let is = i.stone();
                        if is.is_none() {
                            return true
                        } else {
                            if s.is_friend(&is) {
                                i.neighbors().iter().for_each(|n| to_visit.push((n.row, n.col)));
                            }
                        }
                    }
                }
                false // if the search didn't find any empty intersections, this stone is dead
            },
        }
    }
}