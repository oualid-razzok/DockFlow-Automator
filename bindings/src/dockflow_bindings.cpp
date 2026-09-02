// ============================================================================
// dockflow_bindings - C++ accelerator kernels for DockFlow-Automator.
//
// Provides performance-critical operations used by the analyzer and grid
// box modules:
//
//   * parse_pdbqt_atoms   - fast PDBQT ATOM/HETATM line parsing
//   * grid_box            - bounding box + padding from an (N,3) array
//   * pairwise_min_dist   - per-atom minimum distance between two clouds
//   * min_contacts        - all pairs within a cutoff (contact detection)
//   * direct_rmsd         - RMSD without superposition
//   * kabsch_rmsd         - RMSD after optimal rotation (Horn's quaternion
//                            method with a symmetric 4x4 Jacobi eigensolver)
//   * ligand_efficiency   - affinity per heavy atom
//
// The Python fallbacks in dockflow_core keep results bit-for-bit compatible,
// so this module is strictly an accelerator (install via pip wheel, never
// required for correctness).
// ============================================================================

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

namespace py = pybind11;

using ArrayD = py::array_t<double, py::array::c_style | py::array::forcecast>;

static const char kVersion[] = "0.1.0";

// ---------------------------------------------------------------------------
// PDBQT parsing
// ---------------------------------------------------------------------------
struct PdbqtAtom {
    int model = 1;
    int serial = 0;
    std::string name;
    std::string resname;
    std::string chain;
    int resseq = 0;
    double x = 0.0, y = 0.0, z = 0.0;
    double charge = 0.0;
    std::string atom_type;
    std::string record_type;
};

static inline double ColumnFloat(const std::string& line, size_t begin, size_t end) {
    if (line.size() <= begin) return 0.0;
    std::string field = line.substr(begin, end - begin);
    try {
        return std::stod(field);
    } catch (...) {
        return 0.0;
    }
}

static inline int ColumnInt(const std::string& line, size_t begin, size_t end) {
    if (line.size() <= begin) return 0;
    std::string field = line.substr(begin, end - begin);
    try {
        return std::stoi(field);
    } catch (...) {
        return 0;
    }
}

static inline std::string ColumnStr(const std::string& line, size_t begin, size_t end) {
    if (line.size() <= begin) return "";
    std::string field = line.substr(begin, end - begin);
    // trim
    const char* ws = " \t\r\n";
    size_t start = field.find_first_not_of(ws);
    if (start == std::string::npos) return "";
    size_t stop = field.find_last_not_of(ws);
    return field.substr(start, stop - start + 1);
}

static std::vector<py::dict> ParsePdbqtAtoms(const std::string& text) {
    std::vector<py::dict> atoms;
    std::istringstream stream(text);
    std::string line;
    int model = 1;
    while (std::getline(stream, line)) {
        if (line.size() < 6) continue;
        const std::string record = line.substr(0, 6);
        if (record.rfind("MODEL", 0) == 0) {
            model = ColumnInt(line, 10, 14);
            if (model == 0) model = 1;
            continue;
        }
        if (record.rfind("ENDMDL", 0) == 0) {
            ++model;
            continue;
        }
        if (record != "ATOM  " && record != "HETATM") continue;
        py::gil_scoped_release release;  // pure C++ section below
        PdbqtAtom atom;
        atom.model = model;
        atom.record_type = record.substr(0, 6);
        atom.serial = ColumnInt(line, 6, 11);
        atom.name = ColumnStr(line, 12, 16);
        atom.resname = ColumnStr(line, 17, 20);
        atom.chain = ColumnStr(line, 21, 22);
        atom.resseq = ColumnInt(line, 22, 26);
        atom.x = ColumnFloat(line, 30, 38);
        atom.y = ColumnFloat(line, 38, 46);
        atom.z = ColumnFloat(line, 46, 54);
        // PDBQT extension: charge at 70-76, AD4 type at 77-79.  Mirror the
        // Python parser exactly, including the element-column fallback for
        // lines shorter than 79 characters.
        if (line.size() >= 79) {
            atom.charge = ColumnFloat(line, 70, 76);
            atom.atom_type = ColumnStr(line, 77, 79);
        } else if (line.size() >= 76) {
            atom.charge = ColumnFloat(line, 70, 76);
        }
        if (atom.atom_type.empty()) {
            atom.atom_type = ColumnStr(line, 76, 78);
            std::transform(atom.atom_type.begin(), atom.atom_type.end(),
                           atom.atom_type.begin(), ::toupper);
        }
        py::gil_scoped_acquire acquire;
        py::dict dict;
        dict["model"] = atom.model;
        dict["serial"] = atom.serial;
        dict["name"] = atom.name;
        dict["resname"] = atom.resname;
        dict["chain"] = atom.chain;
        dict["resseq"] = atom.resseq;
        dict["x"] = atom.x;
        dict["y"] = atom.y;
        dict["z"] = atom.z;
        dict["charge"] = atom.charge;
        dict["atom_type"] = atom.atom_type;
        dict["record_type"] = atom.record_type.substr(0, record.find_first_of(' '));
        if (dict["record_type"].cast<std::string>().empty()) {
            dict["record_type"] = "ATOM";
        }
        atoms.push_back(std::move(dict));
    }
    return atoms;
}

// ---------------------------------------------------------------------------
// Grid box
// ---------------------------------------------------------------------------
static py::dict GridBox(const ArrayD& coords, double padding, double min_size) {
    py::buffer_info info = coords.request();
    if (info.ndim != 2 || info.shape[1] != 3) {
        throw std::invalid_argument("coords must be an (N, 3) array");
    }
    const size_t n = static_cast<size_t>(info.shape[0]);
    if (n == 0) {
        throw std::invalid_argument("coords is empty");
    }
    const double* data = static_cast<const double*>(info.ptr);
    double lo[3] = {data[0], data[1], data[2]};
    double hi[3] = {data[0], data[1], data[2]};
    for (size_t i = 1; i < n; ++i) {
        for (int d = 0; d < 3; ++d) {
            const double value = data[i * 3 + d];
            if (value < lo[d]) lo[d] = value;
            if (value > hi[d]) hi[d] = value;
        }
    }
    double center[3], size[3];
    for (int d = 0; d < 3; ++d) {
        center[d] = 0.5 * (lo[d] + hi[d]);
        size[d] = std::max(hi[d] - lo[d] + 2.0 * padding, min_size);
    }
    py::dict result;
    result["center"] = std::vector<double>{center[0], center[1], center[2]};
    result["size"] = std::vector<double>{size[0], size[1], size[2]};
    return result;
}

// ---------------------------------------------------------------------------
// Distance kernels
// ---------------------------------------------------------------------------
static ArrayD PairwiseMinDist(const ArrayD& a, const ArrayD& b) {
    py::buffer_info fa = a.request();
    py::buffer_info fb = b.request();
    if (fa.ndim != 2 || fa.shape[1] != 3 || fb.ndim != 2 || fb.shape[1] != 3) {
        throw std::invalid_argument("inputs must be (N, 3) and (M, 3) arrays");
    }
    const size_t na = static_cast<size_t>(fa.shape[0]);
    const size_t nb = static_cast<size_t>(fb.shape[0]);
    if (na == 0 || nb == 0) {
        throw std::invalid_argument("empty coordinate array");
    }
    const double* pa = static_cast<const double*>(fa.ptr);
    const double* pb = static_cast<const double*>(fb.ptr);
    ArrayD output(na);
    py::buffer_info fo = output.request();
    double* po = static_cast<double*>(fo.ptr);
    py::gil_scoped_release release;
    for (size_t i = 0; i < na; ++i) {
        double best = 1e30;
        const double ax = pa[i * 3], ay = pa[i * 3 + 1], az = pa[i * 3 + 2];
        for (size_t j = 0; j < nb; ++j) {
            const double dx = ax - pb[j * 3];
            const double dy = ay - pb[j * 3 + 1];
            const double dz = az - pb[j * 3 + 2];
            const double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < best) best = d2;
        }
        po[i] = std::sqrt(best);
    }
    return output;
}

static std::vector<std::tuple<int, int, double>> MinContacts(
        const ArrayD& a, const ArrayD& b, double cutoff) {
    py::buffer_info fa = a.request();
    py::buffer_info fb = b.request();
    if (fa.ndim != 2 || fa.shape[1] != 3 || fb.ndim != 2 || fb.shape[1] != 3) {
        throw std::invalid_argument("inputs must be (N, 3) and (M, 3) arrays");
    }
    const size_t na = static_cast<size_t>(fa.shape[0]);
    const size_t nb = static_cast<size_t>(fb.shape[0]);
    const double* pa = static_cast<const double*>(fa.ptr);
    const double* pb = static_cast<const double*>(fb.ptr);
    const double cutoff2 = cutoff * cutoff;
    std::vector<std::tuple<int, int, double>> contacts;
    {
        py::gil_scoped_release release;
        for (size_t i = 0; i < na; ++i) {
            const double ax = pa[i * 3], ay = pa[i * 3 + 1], az = pa[i * 3 + 2];
            for (size_t j = 0; j < nb; ++j) {
                const double dx = ax - pb[j * 3];
                const double dy = ay - pb[j * 3 + 1];
                const double dz = az - pb[j * 3 + 2];
                const double d2 = dx * dx + dy * dy + dz * dz;
                if (d2 <= cutoff2) {
                    contacts.emplace_back(static_cast<int>(i), static_cast<int>(j),
                                          std::sqrt(d2));
                }
            }
        }
    }
    std::sort(contacts.begin(), contacts.end(),
              [](const std::tuple<int, int, double>& lhs,
                 const std::tuple<int, int, double>& rhs) {
                  if (std::get<0>(lhs) != std::get<0>(rhs)) {
                      return std::get<0>(lhs) < std::get<0>(rhs);
                  }
                  return std::get<2>(lhs) < std::get<2>(rhs);
              });
    return contacts;
}

// ---------------------------------------------------------------------------
// RMSD
// ---------------------------------------------------------------------------
static double DirectRmsd(const ArrayD& a, const ArrayD& b) {
    py::buffer_info fa = a.request();
    py::buffer_info fb = b.request();
    if (fa.ndim != 2 || fa.shape[1] != 3 || fb.ndim != 2 || fb.shape[1] != 3 ||
        fa.shape[0] != fb.shape[0]) {
        throw std::invalid_argument("inputs must be two (N, 3) arrays of equal length");
    }
    const size_t n = static_cast<size_t>(fa.shape[0]);
    const double* pa = static_cast<const double*>(fa.ptr);
    const double* pb = static_cast<const double*>(fb.ptr);
    py::gil_scoped_release release;
    double total = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double dx = pa[i * 3] - pb[i * 3];
        const double dy = pa[i * 3 + 1] - pb[i * 3 + 1];
        const double dz = pa[i * 3 + 2] - pb[i * 3 + 2];
        total += dx * dx + dy * dy + dz * dz;
    }
    return std::sqrt(total / static_cast<double>(n));
}

// --- Horn's quaternion method -----------------------------------------------
// Builds the symmetric 4x4 N matrix from the cross-covariance and solves the
// eigenproblem with cyclic Jacobi rotations (no external linear algebra
// dependency needed).
static void CenterCloud(std::vector<double>& c) {
    const size_t n = c.size() / 3;
    if (n == 0) return;
    double mx = 0, my = 0, mz = 0;
    for (size_t i = 0; i < n; ++i) {
        mx += c[i * 3]; my += c[i * 3 + 1]; mz += c[i * 3 + 2];
    }
    mx /= n; my /= n; mz /= n;
    for (size_t i = 0; i < n; ++i) {
        c[i * 3] -= mx; c[i * 3 + 1] -= my; c[i * 3 + 2] -= mz;
    }
}

static void JacobiEigen4(std::vector<double>& a, std::vector<double>& eigenvalues) {
    // a: 16 values, row-major symmetric 4x4 (destroyed). eigenvalues: 4.
    for (int sweep = 0; sweep < 64; ++sweep) {
        double off = 0.0;
        for (int p = 0; p < 4; ++p) {
            for (int q = p + 1; q < 4; ++q) off += a[p * 4 + q] * a[p * 4 + q];
        }
        if (off < 1e-24) break;
        for (int p = 0; p < 4; ++p) {
            for (int q = p + 1; q < 4; ++q) {
                if (std::fabs(a[p * 4 + q]) < 1e-18) continue;
                const double theta = 0.5 * (a[q * 4 + q] - a[p * 4 + p]) / a[p * 4 + q];
                const double t = (theta >= 0.0 ? 1.0 : -1.0) /
                                 (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0);
                const double s = t * c;
                for (int k = 0; k < 4; ++k) {
                    const double akp = a[k * 4 + p];
                    const double akq = a[k * 4 + q];
                    a[k * 4 + p] = c * akp - s * akq;
                    a[k * 4 + q] = s * akp + c * akq;
                }
                for (int k = 0; k < 4; ++k) {
                    const double apk = a[p * 4 + k];
                    const double aqk = a[q * 4 + k];
                    a[p * 4 + k] = c * apk - s * aqk;
                    a[q * 4 + k] = s * apk + c * aqk;
                }
            }
        }
    }
    for (int i = 0; i < 4; ++i) eigenvalues[i] = a[i * 4 + i];
}

static double KabschRmsd(const ArrayD& a, const ArrayD& b) {
    py::buffer_info fa = a.request();
    py::buffer_info fb = b.request();
    if (fa.ndim != 2 || fa.shape[1] != 3 || fb.ndim != 2 || fb.shape[1] != 3 ||
        fa.shape[0] != fb.shape[0]) {
        throw std::invalid_argument("inputs must be two (N, 3) arrays of equal length");
    }
    const size_t n = static_cast<size_t>(fa.shape[0]);
    if (n < 3) {
        return DirectRmsd(a, b);
    }
    std::vector<double> pa(static_cast<const double*>(fa.ptr),
                           static_cast<const double*>(fa.ptr) + n * 3);
    std::vector<double> pb(static_cast<const double*>(fb.ptr),
                           static_cast<const double*>(fb.ptr) + n * 3);
    py::gil_scoped_release release;
    CenterCloud(pa);
    CenterCloud(pb);
    // Cross-covariance S = sum p_i q_i^T
    double s[3][3] = {{0}};
    for (size_t i = 0; i < n; ++i) {
        for (int d1 = 0; d1 < 3; ++d1) {
            for (int d2 = 0; d2 < 3; ++d2) {
                s[d1][d2] += pa[i * 3 + d1] * pb[i * 3 + d2];
            }
        }
    }
    // Horn's symmetric 4x4 matrix (quaternion parametrisation).
    const double trace = s[0][0] + s[1][1] + s[2][2];
    std::vector<double> horn(16, 0.0);
    horn[0] = trace;
    horn[1] = s[1][2] - s[2][1];
    horn[2] = s[2][0] - s[0][2];
    horn[3] = s[0][1] - s[1][0];
    horn[4] = horn[1];
    horn[5] = 2.0 * s[0][0] - trace;
    horn[6] = s[0][1] + s[1][0];
    horn[7] = s[0][2] + s[2][0];
    horn[8] = horn[2];
    horn[9] = horn[6];
    horn[10] = 2.0 * s[1][1] - trace;
    horn[11] = s[1][2] + s[2][1];
    horn[12] = horn[3];
    horn[13] = horn[7];
    horn[14] = horn[11];
    horn[15] = 2.0 * s[2][2] - trace;
    std::vector<double> eigenvalues(4, 0.0);
    JacobiEigen4(horn, eigenvalues);
    const double max_eigenvalue = *std::max_element(eigenvalues.begin(), eigenvalues.end());
    // rmsd^2 = (|P|^2 + |Q|^2 - 2 lambda_max) / N
    double norm = 0.0;
    for (size_t i = 0; i < n * 3; ++i) {
        norm += pa[i] * pa[i] + pb[i] * pb[i];
    }
    const double value = (norm - 2.0 * max_eigenvalue) / static_cast<double>(n);
    return std::sqrt(std::max(value, 0.0));
}

// ---------------------------------------------------------------------------
// Misc
// ---------------------------------------------------------------------------
static double LigandEfficiency(double affinity, int num_heavy_atoms) {
    if (num_heavy_atoms <= 0) {
        throw std::invalid_argument("num_heavy_atoms must be positive");
    }
    return affinity / static_cast<double>(num_heavy_atoms);
}

static std::vector<std::vector<double>> BoxCorners(const ArrayD& center,
                                                   const ArrayD& size) {
    py::buffer_info fc = center.request();
    py::buffer_info fs = size.request();
    if (fc.ndim != 1 || fc.shape[0] != 3 || fs.ndim != 1 || fs.shape[0] != 3) {
        throw std::invalid_argument("center and size must be 3-vectors");
    }
    const double* pc = static_cast<const double*>(fc.ptr);
    const double* ps = static_cast<const double*>(fs.ptr);
    std::vector<std::vector<double>> corners(8, std::vector<double>(3));
    int index = 0;
    for (int sx = 0; sx < 2; ++sx) {
        for (int sy = 0; sy < 2; ++sy) {
            for (int sz = 0; sz < 2; ++sz) {
                corners[index][0] = pc[0] + (sx ? 0.5 : -0.5) * ps[0];
                corners[index][1] = pc[1] + (sy ? 0.5 : -0.5) * ps[1];
                corners[index][2] = pc[2] + (sz ? 0.5 : -0.5) * ps[2];
                ++index;
            }
        }
    }
    return corners;
}

// ---------------------------------------------------------------------------
// Module
// ---------------------------------------------------------------------------
PYBIND11_MODULE(dockflow_bindings, m) {
    m.doc() = "DockFlow-Automator C++ accelerator kernels (PDBQT parsing, "
              "grid box, contacts, RMSD).";
    m.attr("__version__") = kVersion;

    m.def("parse_pdbqt_atoms", &ParsePdbqtAtoms, py::arg("text"),
          "Parse PDBQT ATOM/HETATM records into a list of dicts.");

    m.def("grid_box", &GridBox, py::arg("coords"), py::arg("padding") = 4.0,
          py::arg("min_size") = 10.0,
          "Bounding box + padding from an (N,3) coordinate array.");

    m.def("pairwise_min_dist", &PairwiseMinDist, py::arg("a"), py::arg("b"),
          "Per-row minimum distance from cloud a to cloud b.");

    m.def("min_contacts", &MinContacts, py::arg("a"), py::arg("b"),
          py::arg("cutoff") = 5.0,
          "All (i, j, distance) pairs between two clouds within a cutoff.");

    m.def("direct_rmsd", &DirectRmsd, py::arg("a"), py::arg("b"),
          "RMSD without superposition.");

    m.def("kabsch_rmsd", &KabschRmsd, py::arg("a"), py::arg("b"),
          "RMSD after optimal rigid superposition (Horn's quaternion method).");

    m.def("ligand_efficiency", &LigandEfficiency, py::arg("affinity"),
          py::arg("num_heavy_atoms"), "Affinity per heavy atom.");

    m.def("box_corners", &BoxCorners, py::arg("center"), py::arg("size"),
          "The 8 corners of a grid box.");
}
