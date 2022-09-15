using Random, LinearAlgebra, Distributions, StatsBase
# Linear Algebra of Markov chain

# Example: Weather Markov chain of Oz
S = [:R, :S, :C]
P = [0.5 0.25 0.25
     0.5  0    0.5
     0.25 0.25 0.5]

sum(P[2, :])

display(P*P)

P2 = [0.4375  0.1875  0.375
      0.375   0.25    0.375
      0.375   0.1875  0.4375]

# Probability vector for today
p = [1/3 2/3 0]

# Probability vector for the future
p2 = p * P#[0.5 0.0833333 0.41]
p3 = p2 * P

@show p
@show p2
@show p3

# Associativity of matrix multiplication
p
p2 = p*P
p3 = p2*P
p3 ≈ p * (P*P) # p3 and p*(P*P) are numerically the same

Q = [0.4 0.2 0.4
     0.4 0.2 0.4
     0.4 0.2 0.4]
q = [0.4 0.2 0.4]

# start i = 2
# what is the chance that i am in state 3 after 10 steps

(P^10)[2, 3] ≈ Q[2, 3]
Q[2,3] ≈ Q[1,3] ≈ Q[3,3]
q[3]

# Compute q 
qvector(P) = let x = nullspace((I-P)'); x/sum(x); end # works, but not most elegant way of computing eigenvectors...
qvector(P)


# Implementation: how to sample a Markov chain
"""
    samplefrom(p)

Produces a sample from the probability distribution
with prob .mass vector p.

Helpful: to sample from the conditional distribution given state i

    j = samplefrom(P[i, :]) # : for all j

"""
samplefrom(p) = sample(1:3, weights(p))

samplefrom(p)
p

# starting vector
Random.seed!(1)
s0 = samplefrom(p)

s1 = samplefrom(P[s0, :])

s2 = samplefrom(P[s1, :])

s3 = samplefrom(P[s2, :])

[s0, s1, s2, s3]

"""
    samplefromchain(p, P, n)

Samples an n-step Markov chain with starting distribution p
with transition matrix P.
"""
function samplefromchain(p, P, n) # n is number of trans
    s = samplefrom(p)
    statevector = [s]
    for i in 1:n
        s = samplefrom(P[s, :])
        push!(statevector, s)
    end
    return statevector
end

chain = samplefromchain(p, P, 3)

chain

map(i -> ["R", "S", "C"][i], chain) # map over the keys


# example absorbing:

# Transition matrix


        #A  #B #C #D
P = [   1   0   0   0       # A absorbing state
        0   .5  0.5 0       # B
        0   0   1   0       # C absorbing state
        1/4 1/2 0 1/4]      # D

# 3 is absorbing because
1 - P[3,3]
#the probability to leave 3 is zero

i = 2
check_absorbing(P, i) = (P[i, i] == 1)


# Order the elements
absorbing = [check_absorbing(P, i) for i in 1:size(P,1)]

order = sortperm(absorbing)
println("New order:", ('A':'D')[order])
Px = [P[i,j] for i in order, j in order]
display(Px)
