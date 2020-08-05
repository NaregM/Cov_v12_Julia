using Pkg
using LinearAlgebra
using Distributions
using FFTW
using DelimitedFiles
using PyPlot
using Plots
using Random
using LaTeXStrings
using CSV, DataFrames
using ProgressMeter

println("Imported Libs: OK")

# Functions ============================
function v12(df, M_cut)

    M200c = df[:, 1];
    ind = @. M200c * 1 /1 > M_cut # 1 -> 1e10 for E

    M = M200c[ind]

    nclus = length(M200c[ind])
    println(nclus)

    pos = df[ind, 2:4]./1.0/1   # /1 -> /1000
    #pos[:, 1] .+= 1000.0
    vel = df[ind, 5:7]

    nbins = 40
    binsize = 4.0
    rbins = binsize * (range(0, stop = nbins - 1, step = 1) .+ 0.5)

    v_of_r = zeros(nbins)
    n_of_r = zeros(nbins)

    v12_full = zeros(Float32, nclus, nclus)
    ibin_full = zeros(Int32,(nclus, nclus))

    for i in 1:nclus, j in 1:i

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

            v12_full[i, j] = dot(t1 .- t2, q_AB)#dv_pair_subf

            if ibin <= nbins

                ibin_full[i, j] = ibin

            end

        end

    end

    return rbins, v12_full, ibin_full

end


function v12_perbin(df, n)

    if n in df[:, 2]


        return mean(df[:, 1][df[:, 2] .== n])

    else

        return 0.0

    end

end


function divide_sample(df, n_sample)

    N = size(df)[1]

    subsample = [df[i: min(i + n_sample - 1, end), :] for i in 1:n_sample:N];

    data_jk = []

    for k in 1:length(subsample)

        append!(data_jk, [vcat([subsample[i] for i = 1:length(subsample) if i != k]...)])

    end

    return data_jk

end
# ======================================

df = readdlm("/media/nareg/My Passport/Quijote_Sim/Catalogs/cat1.txt")
df1 = readdlm("/home/nareg/Downloads/websky_BG_SOxDESI.txt")
cat_i = "_websky"


N_jk = 20;
M_thr = 7.0e14;

r_, v12FULL_, binFULL_ = v12(df, M_thr);

# Preparing data
n_r = length(r_);

binFULL = binFULL_[binFULL_ .!= 0.0];
v12FULL = v12FULL_[binFULL_ .!= 0.0];

N_full = length(binFULL);

v12_bin = hcat(v12FULL, binFULL)   #[[v12FULL[i]; binFULL[i]] for i in 1:N_full]; #
v12_jk = divide_sample(v12_bin, Int16(ceil(N_full/N_jk)))
# ==================================================================

# Covariance Matrix
C = zeros(n_r, n_r)

@showprogress 1 "Computing..." for ri in 2:n_r-1, rj in 2:n_r-1

    C_ = 0


    for n in 1: N_jk

        C_ += (-v12_perbin(v12_jk[n], ri) - -v12_perbin(v12_bin, ri)) * (-v12_perbin(v12_jk[n], rj) - -v12_perbin(v12_bin, rj))

    end

    C[ri, rj] = C_

end



writedlm("Cov_v12$cat_i.txt", C)

println("Finished! Results Saved @ cov_v12$cat_i.txt")
