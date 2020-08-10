using Pkg
using LinearAlgebra
using Distributions
using DelimitedFiles
using PyPlot
using Plots
using Random
using ProgressMeter

# ==========================================

# pairwise velocity among all posibile pairs
function v12(df, M_cut, v12_type)

    M200c = df[:, 1];
    ind = @. M200c * 1 /1 > M_cut # 1 -> 1e10 for E

    M = M200c[ind]

    nclus = length(M200c[ind])
    println(nclus)

    pos = df[ind, 2:4]./1.0/1   # /1 -> /1000
    #pos[:, 1] .+= 1000.0
    vel = df[ind, 5:7]

    nbins = 30
    binsize = 4.0
    rbins = binsize * (range(0, stop = nbins - 1, step = 1) .+ 0.5)

    v_of_r = zeros(nbins)
    n_of_r = zeros(nbins)

    v12_full = zeros(Float32, nclus, nclus)
    ibin_full = zeros(Int32,(nclus, nclus))

    @showprogress "Initializing ..." for i in 1:nclus, j in 1:i

        if (i != j)

            dr = pos[i,:] - pos[j,:]

            r1 = pos[i,:]./norm(pos[i,:])
            r2 = pos[j,:]./norm(pos[j,:])

            vr1 = dot(vel[i,:], r1)
            vr2 = dot(vel[j,:], r2)

            drabs = sqrt(dot(dr, dr))
            dr_norm = dr/drabs

            dv_subf = vel[i,:] - vel[j,:]
            dv_pair_subf = dot(dv_subf, dr_norm)

            ibin = Int(floor(drabs/binsize)) + 1

            q_AB = 0.5 * (2 * dr_norm - dot(dr_norm, r1)*r1 - dot(dr_norm, r2)*r2)
            t1 = vel[i,:]# - vr1*r1
            t2 = vel[j,:]# - vr2*r2

            if v12_type == "transverse"

                v12_full[i, j] = dot(t1 .- t2, q_AB)#dv_pair_subf

            elseif v12_type == "direct"

                v12_full[i, j] = dv_pair_subf

            elseif v12_type == "radial"

                v12_full[i, j] = dv_pair_subf

            else

                println("Wrong Choice for v12 Estimator")
                break

            end


            if ibin <= nbins

                ibin_full[i, j] = ibin

            end

        end

    end

    return rbins, v12_full, ibin_full

end


# Calculate the average pairwise velocity of each seperation bin
# Each element of df is a tuple (a, b) where a is pairwise velocity
# and b is separation-bin it belongs to
function v12_perbin(df, n)

    if n in df[:, 2]


        return mean(df[:, 1][df[:, 2] .== n])

    else

        return 0.0

    end

end


# Number of pairs in seperation bin
function n_perbin(df, n)

    if n in df[:, 2]


        return sum(v12_bin[:, 2] .== n)

    else

        return 0.0

    end

end


# divide df to length(df)/n_sample subsamples
# Take one possible subsample of size (length(df)/n_sample) out at a time and unite the rest
# return all possible combinations
function divide_sample(df, n_sample)

    N = size(df)[1]

    subsample = [df[i: min(i + n_sample - 1, end), :] for i in 1:n_sample:N];

    data_jk = []

    for k in 1:length(subsample)

        append!(data_jk, [vcat([subsample[i] for i = 1:length(subsample) if i != k]...)])

    end

    return data_jk

end
