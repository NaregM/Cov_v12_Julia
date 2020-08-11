include("julia_v12_cov.jl")


#df = readdlm("/media/nareg/My Passport/Quijote_Sim/Catalogs/cat1.txt")
#df2 = readdlm("/media/nareg/My Passport/Quijote_Sim/Catalogs/cat1.txt")
df = readdlm("/home/nareg/Downloads//websky_BG_SOxDESI_txmap.txt")

cat_i = "_websky"

# Number of separation bins, separation bin size and threshold mass
binsize = 4.0
nbin = 50
M_thr = 0.8e15

N_jk = 100;

# r: is the separation between halos
# v12_matrix: paiwise velocity for all possible pairs
# bin_matrix: separation bin of the pair
r, v12_matrix, bin_matrix = v12(df, M_thr, "direct");

# Preparing data
# ====================================================================

n_r = length(r)
# flatten; throw away zero elements
bins_ = bin_matrix[bin_matrix .!= 0.0]
v12s = v12_matrix[bin_matrix .!= 0.0]

# Number of all possible pairs
N = length(bins_)

# Create a list of (pariwise, bin) list
v12_bin = hcat(v12s, bins_)

# Divide the sameple into JK samples
v12_JK = divide_sample(v12_bin, Int16(ceil(N/N_jk)))  # create 100 JK samples

# Average pariwise velocity of each separation-bin
v12_mean = [v12_perbin(v12_bin, ri) for ri in 2:n_r-1]

# ==================================================================

# Covariance Matrix
C = zeros(n_r, n_r)

@showprogress 1 "Computing..." for ri in 2:n_r-1, rj in 2:n_r-1

    C_ = 0


    for n in 1: N_jk

        C_ += (-v12_perbin(v12_JK[n], ri) - -v12_perbin(v12_bin, ri)) * (-v12_perbin(v12_JK[n], rj) - -v12_perbin(v12_bin, rj))

    end

    C[ri, rj] = C_

end

# ======================================================================
# Save the results
writedlm("SNR_$cat_i.txt", SNR(v12_mean, C))
writedlm("Cov_v12_$cat_i.txt", C)
writedlm("v12_avg_$cat_i.txt", v12_mean)

println("Finished! Results Saved @ SNR_$cat_i.txt")
println("Finished! Results Saved @ Cov_v12_$cat_i.txt")
println("Finished! Results Saved @ v12_avg_$cat_i.txt")








