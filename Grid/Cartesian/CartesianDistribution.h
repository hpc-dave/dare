/*
 * MIT License
 *
 * Copyright (c) 2024 David Rieder

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef GRID_CARTESIAN_CARTESIANDISTRIBUTION_H_
#define GRID_CARTESIAN_CARTESIANDISTRIBUTION_H_
#include <algorithm>
#include <cmath>
#include <vector>

#include "MPI/ExecutionManager.h"
#include "Utilities/Vector.h"

namespace dare::Grid {

namespace details {
/*!
 * \brief routine for decomposition for Cartesian grid
 * This decomposition aims at providing as many cubical domains as possible
 * and will only subdivide the domains at the longest end unequally.
 */

/*!
 * \brief routine for decomposition for Cartesian grid
 * This decomposition aims at providing as many cubical domains as possible
 * and will only subdivide the domains at the longest end unequally.
 */
template <std::size_t Dim, class LO, class GO>
void CartesianDistribution_MPI_Dims_create(int num_proc,
                                  const utils::Vector<Dim, GO>& resolution_global,
                                  std::vector<utils::Vector<Dim, LO>>* vec_res_local,
                                  std::vector<utils::Vector<Dim, GO>>* vec_offsets) {
    vec_res_local->resize(num_proc);
    vec_offsets->resize(num_proc);
    int ndims = Dim;
    int dims[Dim];
    std::fill(dims, dims + Dim, 0);
    MPI_Dims_create(num_proc, ndims, dims);
    utils::Vector<Dim, GO> subdomain_res, subdomain_add_to_last;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        subdomain_res[dim] = resolution_global[dim] / dims[dim];
        subdomain_add_to_last[dim] = resolution_global[dim] - subdomain_res[dim] * dims[dim];
    }

    utils::Vector<Dim, int> hsum_topo;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        hsum_topo[dim] = 1;
        for (std::size_t n{dim + 1}; n < Dim; n++)
            hsum_topo[dim] *= dims[n];
    }

    auto GetIndex = [&](int n) {
        utils::Vector<Dim, int> ind;
        for (std::size_t dim{0}; dim < Dim; dim++) {
            ind[dim] = n / hsum_topo[dim];
            n -= ind[dim] * hsum_topo[dim];
        }
        return ind;
    };

    for (int n_proc{0}; n_proc < num_proc; n_proc++) {
        (*vec_res_local)[n_proc] = subdomain_res;
        utils::Vector<Dim, int> ind = GetIndex(n_proc);
        for (std::size_t dim{0}; dim < Dim; dim++) {
            if (ind[dim] == (dims[dim] - 1)) {
                (*vec_res_local)[n_proc][dim] += subdomain_add_to_last[dim];
            }
            (*vec_offsets)[n_proc][dim] = ind[dim] * subdomain_res[dim];
        }
    }
}

template <std::size_t Dim, class LO, class GO>
void CartesianDistribution_Cubical(int num_proc,
                                  const utils::Vector<Dim, GO>& resolution_global,
                                  std::vector<utils::Vector<Dim, LO>>* vec_res_local,
                                  std::vector<utils::Vector<Dim, GO>>* vec_offsets, bool print_warning = true) {
    const double max_ratio_deviation{1.2};  // acceptable volume ration of domains
    // Determine direction with maximum number of cells
    vec_res_local->resize(num_proc);
    vec_offsets->resize(num_proc);
    std::size_t pos_max_res{0};
    for (std::size_t n{1}; n < Dim; n++)
        if (resolution_global[n] > resolution_global[pos_max_res])
            pos_max_res = n;

    // Determine approximately number of cells per subdomain
    std::size_t num_cells_total{1};
    for (auto e : resolution_global)
        num_cells_total *= e;
    double avg_cells_proc = static_cast<double>(num_cells_total) / num_proc;
    double length_side = std::pow(avg_cells_proc, 1. / Dim);

    // set and correct number of cells per direction of subdomain
    dare::utils::Vector<Dim, GO> subdomain_res;
    subdomain_res.SetAllValues(static_cast<GO>(length_side));

    for (std::size_t dim{0}; dim < Dim; dim++)
        subdomain_res[dim] = std::min(subdomain_res[dim], resolution_global[dim]);

    dare::utils::Vector<Dim, int> topo, topo_hierarchic_sum;
    // remaining cells which cannot be evenly distributed
    dare::utils::Vector<Dim, GO> remaining_cells;

    for (std::size_t dim{0}; dim < Dim; dim++) {
        topo[dim] = resolution_global[dim] / subdomain_res[dim];
        remaining_cells[dim] = resolution_global[dim] - (subdomain_res[dim] * topo[dim]);
        subdomain_res[dim] += remaining_cells[dim] / topo[dim];
        if (dim == pos_max_res)
            remaining_cells[dim] = 0;  // special treatment later
        else
            remaining_cells[dim] = resolution_global[dim] - (subdomain_res[dim] * topo[dim]);
    }
    int delta_proc = 1;
    for (auto e : topo)
        delta_proc *= e;
    delta_proc = num_proc - delta_proc;
    int num_proc_slice = 1;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        if (dim == pos_max_res)
            continue;
        num_proc_slice *= topo[dim];
    }
    int missing_slices = delta_proc / num_proc_slice;
    if (delta_proc > 0 && (delta_proc % num_proc_slice))
        missing_slices += 1;

    topo[pos_max_res] += missing_slices;
    subdomain_res[pos_max_res] = resolution_global[pos_max_res] / topo[pos_max_res];
    remaining_cells[pos_max_res] = resolution_global[pos_max_res] - (subdomain_res[pos_max_res] * topo[pos_max_res]);

    for (std::size_t dim{0}; dim < Dim; dim++) {
        topo_hierarchic_sum[dim] = 1;
        for (std::size_t n{dim + 1}; n < Dim; n++)
            topo_hierarchic_sum[dim] *= topo[n];
    }

    for (int n{0}; n < num_proc; n++) {
        utils::Vector<Dim, int> proc_ind;
        int loc_n{n};
        for (std::size_t dim{0}; dim < Dim; dim++) {
            proc_ind[dim] = loc_n / topo_hierarchic_sum[dim];
            loc_n -= proc_ind[dim] * topo_hierarchic_sum[dim];
        }
        for (std::size_t dim{0}; dim < Dim; dim++) {
            (*vec_res_local)[n][dim] = subdomain_res[dim];
            (*vec_offsets)[n][dim] = subdomain_res[dim] * proc_ind[dim];
            if (proc_ind[dim] == (topo[dim] - 1)) {
                (*vec_res_local)[n][dim] += remaining_cells[dim];
            }
        }
    }

    // at this point we have equally distributed domains, with
    // a few missing cells in the longest direction, those need to be distributed now
    // among the last layer
    // To do so, we determine the number of processes in the last slice and then subdivide
    int n_subdomain = topo[0];
    for (std::size_t dim{1}; dim < Dim; dim++)
        n_subdomain *= topo[dim];

    int surplus_domains = n_subdomain - num_proc;

    int n_pre_last_slice{1};
    for (std::size_t dim{0}; dim < Dim; dim++) {
        if (dim == pos_max_res) {
            n_pre_last_slice *= topo[dim] - 1;
        } else {
            n_pre_last_slice *= topo[dim];
        }
    }

    // // Note, that in 1D, no remaining splitting can occur, so we don't need to care for that case
    int mod{0}, counter{0};
    std::size_t slice_dir{0};
    if (slice_dir == pos_max_res)
        slice_dir++;

    // Here we need to reduce the number of subdomains as indicated by surplus domains
    // For this purpose, all subdomains in the last slice in the dominant direction
    // are merged to 1 and subsequently split into two, until we have created all required
    // subdomains
    bool failed{false};  // intermediate feature, until the issues below are resolved
    for (int n{0}; n < surplus_domains; n++) {
        // for the first processor in the final slice, we reset it to account for the whole slice
        if (n == 0) {
            for (std::size_t dim{0}; dim < Dim; dim++) {
                if (dim != pos_max_res) {
                    (*vec_res_local)[n_pre_last_slice + n][dim] = resolution_global[dim];
                    (*vec_offsets)[n_pre_last_slice + n][dim] = 0;
                }
            }
            mod = 2;
        } else {
            // Subsequently, all elements in the slice are halfed
            // Once we have done that with all the elements in the slice
            // we start again, until we have reached the number of
            // desired processes
            if (static_cast<std::size_t>(n_pre_last_slice + counter) >= vec_res_local->size()) {
                failed = true;
                break;
            }
            (*vec_res_local)[n] = (*vec_res_local)[n_pre_last_slice + counter];
            (*vec_offsets)[n] = (*vec_offsets)[n_pre_last_slice + counter];

            GO subres = (*vec_res_local)[n_pre_last_slice + counter][slice_dir];
            GO subres_first = subres / 2;
            GO subres_second = subres - subres_first;
            (*vec_res_local)[n_pre_last_slice + counter][slice_dir] = subres_first;
            (*vec_res_local)[n][slice_dir] = subres_second;
            (*vec_offsets)[n][slice_dir] = (*vec_offsets)[n_pre_last_slice + counter][slice_dir] + subres_first;

            if (counter == mod) {
                counter = 0;
                mod *= 2;
                slice_dir++;
                if (slice_dir == pos_max_res)
                    slice_dir++;
                if (slice_dir == Dim) {
                    slice_dir = 0;
                    if (pos_max_res == 0)
                        slice_dir = 1;
                }
            } else {
                counter++;
            }
        }
    }

    if (!failed) {
        // control the max volume occupied by subdomains on the last slice
        GO num_cells_regular{1};
        for (auto e : (*vec_res_local)[0])
            num_cells_regular *= e;

        double max_ratio{1.};
        GO delta_cell_max{0};
        for (int n{n_pre_last_slice}; n < num_proc; n++) {
            GO num_cells_sd{1};
            for (auto e : (*vec_res_local)[n])
                num_cells_sd *= e;
            max_ratio = std::max(max_ratio, static_cast<double>(num_cells_sd) / num_cells_regular);
            delta_cell_max = std::max(num_cells_sd - num_cells_regular, delta_cell_max);
        }

        // if the inequality is too big, we remove some layers of the last slice and distribute them
        // on the previous ones
        if ((topo[pos_max_res] > 1) && (max_ratio > max_ratio_deviation)) {
            GO cell_per_slice{1};
            for (std::size_t dim{0}; dim < Dim; dim++)
                if (dim != pos_max_res)
                    cell_per_slice *= resolution_global[dim];

            // total number of cells to redistribute
            LO n_cell_redistribute = delta_cell_max;
            // number of cells in main direction, which will be removed per subdomain slice
            LO n_cell_redist_slice = std::max(GO(1), delta_cell_max / cell_per_slice / (topo[pos_max_res] - 1));

            int proc_per_slice{1};
            for (std::size_t dim{0}; dim < Dim; dim++) {
                if (dim == pos_max_res)
                    continue;
                proc_per_slice *= topo[dim];
            }
            int slice_pre_last_slice = topo[pos_max_res] - 2;
            for (int n{slice_pre_last_slice}; n >= 0; n--) {
                if (n_cell_redistribute <= 0)
                    break;
                int proc_start = n * proc_per_slice;

                for (int n_proc{proc_start}; n_proc < (proc_start + proc_per_slice); n_proc++)
                    (*vec_res_local)[n_proc][pos_max_res] += n_cell_redist_slice;

                for (int n_proc{proc_start + proc_per_slice}; n_proc < n_pre_last_slice; n_proc++)
                    (*vec_offsets)[n_proc][pos_max_res] += n_cell_redist_slice;

                for (int n_proc{n_pre_last_slice}; n_proc < num_proc; n_proc++) {
                    (*vec_offsets)[n_proc][pos_max_res] += n_cell_redist_slice;
                    (*vec_res_local)[n_proc][pos_max_res] -= n_cell_redist_slice;
                }
                n_cell_redistribute -= n_cell_redist_slice * cell_per_slice;
            }
        }
    }

    // if a bug is detected, switch to less efficient MPI_Dims_create
    if (!failed) {
        num_cells_total = 1;
        for (auto e : resolution_global)
            num_cells_total *= e;
        std::size_t sum_cells{0};
        for (const auto& r_sub : *vec_res_local) {
            LO n_sub{1};
            for (LO dim : r_sub)
                n_sub *= dim;
            sum_cells += n_sub;
        }
        failed = sum_cells != num_cells_total;
    }
    if (failed) {
        if (print_warning)
            std::cout << "Cubical Cartesian distribution failed, switching to MPI_Dims_create instead!" << std::endl;
        details::CartesianDistribution_MPI_Dims_create(num_proc, resolution_global, vec_res_local, vec_offsets);
    }
}

}  // namespace details

/*!
 * \brief decomposition for Cartesian grid
 * This decomposition aims at providing as many cubical domains as possible
 * and will only subdivide the domains at the longest end unequally.
 */
template <std::size_t Dim, class LO, class GO>
void CartesianDistribution_Cubical(mpi::ExecutionManager* exec_man,
                                  const utils::Vector<Dim, GO>& resolution_global,
                                  utils::Vector<Dim, LO>* resolution_local,
                                  utils::Vector<Dim, GO>* offset) {
    int tag_res{1000};                      // tag for MPI communication
    int tag_off{1001};                      // tag for MPI communication

    if (exec_man->AmIRoot()) {
        int num_proc = exec_man->GetNumberProcesses();

        // storage for each subdomain
        std::vector<utils::Vector<Dim, LO>> vec_res_local(num_proc);
        std::vector<utils::Vector<Dim, GO>> vec_offsets(num_proc);

        details::CartesianDistribution_Cubical(num_proc, resolution_global,
                                              &vec_res_local, &vec_offsets);

        if (!exec_man->IsSerial()) {
            // communicate the result with remaining processes
            std::vector<MPI_Request> requests_res(num_proc - 1), requests_off(num_proc - 1);
            int count{0};
            for (int n{0}; n < num_proc; n++) {
                if (n == exec_man->GetRank())
                    continue;
                exec_man->Isend(vec_res_local[n].data(), vec_res_local[n].size(), n, tag_res, &requests_res[count]);
                exec_man->Isend(vec_offsets[n].data(), vec_offsets[n].size(), n, tag_off, &requests_off[count]);
                count++;
            }

            std::vector<MPI_Status> status_res(num_proc - 1), status_off(num_proc - 1);
            int status_res_all = MPI_Waitall(requests_res.size(), requests_res.data(), status_res.data());
            int status_off_all = MPI_Waitall(requests_off.size(), requests_off.data(), status_off.data());
            if (status_res_all != MPI_SUCCESS) {
                exec_man->Print(dare::mpi::Verbosity::Low)
                    << "An error occured during the sending of resolutions: " << status_res_all << std::endl;
            }
            if (status_off_all != MPI_SUCCESS) {
                exec_man->Print(dare::mpi::Verbosity::Low)
                    << "An error occured during the sending of offsets: " << status_off_all << std::endl;
            }
        }
        *resolution_local = vec_res_local[exec_man->GetRank()];
        *offset = vec_offsets[exec_man->GetRank()];
    } else {
        exec_man->Recv(resolution_local->data(), resolution_local->size(), exec_man->GetRankRoot(), tag_res);
        exec_man->Recv(offset->data(), offset->size(), exec_man->GetRankRoot(), tag_off);
    }
}

/*!
 * \brief decomposition for Cartesian grid
 * Employs MPI_Dims_create
 */
template <std::size_t Dim, class LO, class GO>
void CartesianDistribution_MPI_Dims_create(mpi::ExecutionManager* exec_man,
                                  const utils::Vector<Dim, GO>& resolution_global,
                                  utils::Vector<Dim, LO>* resolution_local,
                                  utils::Vector<Dim, GO>* offset) {
    int num_proc = exec_man->GetNumberProcesses();

    // storage for each subdomain
    std::vector<utils::Vector<Dim, LO>> vec_res_local(num_proc);
    std::vector<utils::Vector<Dim, GO>> vec_offsets(num_proc);

    details::CartesianDistribution_MPI_Dims_create(num_proc, resolution_global,
                                                   &vec_res_local, &vec_offsets);

    *resolution_local = vec_res_local[exec_man->GetRank()];
    *offset = vec_offsets[exec_man->GetRank()];
}

}  // namespace dare::Grid

#endif  // GRID_CARTESIAN_CARTESIANDISTRIBUTION_H_
