#using
using LinearAlgebra
using QuantEcon

#returns probability of a particular instance (j initial, i end) after n steps given a transition matrix
function P(t_matrix, i, j , n)
    X = t_matrix^n
    return X[i,j]
end
A = [.5 .2 .3; .6 .2 .2; .1 .8 .1]
i = 1
j = 3
n = 2
P(A, i, j, n)
#0.22

#most likely walk with n steps
function walk(t_matrix, n)
    max = maximum(t_matrix)
    first_step = findall(x->x==max,t_matrix)[1]
    index = first_step[2]
    answer = []
    if n == 1
        answer = [first_step[1],first_step[2]]
        return answer
    end
    answer = [first_step[1],first_step[2]]
    if n > 1
        for i = 1:n-1
            p_matrix = t_matrix^i
            max = maximum(p_matrix[index])
            index = (findall(x->x==max,p_matrix)[1])[2]
            answer = push!(answer, index)
        end
    end
    print(join(answer, "->"))
end 
A = [.5 .2 .3; .6 .2 .2; .1 .8 .1]
n = 2
walk(A,10)
#3->2->1->1->1->1->1->1->1->1->1

#returns stationary distribution given transition matrix through ITERATION (n times)
function i_stationary(t_matrix, n)
    x = size(t_matrix)[1]
    station = zeros(x,1)
    for i = 1:x
        initial = zeros(x,1)
        initial[i,1] = 1
        for i = 1:n
            initial = t_matrix*initial
        end
        station[i,1] = initial[1]
    end
    station = reshape(station,1,x)
return station
end
B = [.2 .6 .2 ; .3 0 .7 ; .5 0 .5 ]
i_stationary(B,100)
#=1×3 Matrix{Float64}:
0.352113  0.211268  0.43662=# 

#=returns stationary distribution given transition matrix through eigendecomposition
using QuantEcon=#
function s_stationary(t_matrix)
    mc = MarkovChain(t_matrix)
    answer = stationary_distributions(mc)
    #answer = split(answer)
    return answer
end
B = [.2 .6 .2 ; .3 0 .7 ; .5 0 .5 ]
s_stationary(B)  
#=1-element Vector{Vector{Float64}}:
 [0.35211267605633806, 0.21126760563380287, 0.4366197183098592]=#
    
#returns stationary distribution given transition matrix through eigendecomposition
function d_stationary(t_matrix)
    F = eigen(t_matrix)
    QL = inv(eigvecs(F))
    x = eigvals(F)
    y = size(x)[1]
    ev = 0
    for i = 1:y
        if floor(Int64,real.(x[i])) == 1
            ev = i 
        end
    end
    answer = real.(reshape(QL[ev,1:ev]/sum(QL[ev,1:ev]) ,1,ev))
    return answer
end
B = [.2 .6 .2; .3 0 .7; .5 0 .5]
d_stationary(B)
#= 1×3 Matrix{Float64}:
 0.352113  0.211268  0.43662 =#
    
 #returns probability of observed sequence
function forward_probability(pi, observation, t_matrix, e_matrix)
    N = length(pi)
    T = length(observation)
    forward_matrix = zeros(N,T)
    first_obs = observation[1]
    #initialization
    for i = 1:N
        forward_matrix[i,1] = pi[1,i]*e_matrix[first_obs,i]
    end
    #recursion
    for i = 2:T
        for k = 1:N
            previous = forward_matrix[i-1:i]
            x = previous[1]*t_matrix[1,k]*e_matrix[observation[i],k]+previous[2]*t_matrix[2,k]*e_matrix[observation[i],k]
            forward_matrix[k,i] = x
        end
    end
    #termination
    return sum(forward_matrix[1:N,T])
end
pi = [.8 .2]
e_matrix = [.2 .5; .4 .4; .4 .1]
t_matrix = [.6 .4;.5 .5]
observation = [3 1 3]
forward_probability(pi, observation, t_matrix, e_matrix)
#0.015700000000000002
    
#returns most likely sequence from observed values
function viterbi(pi, observation, t_matrix, e_matrix)
    N = length(pi)
    T = length(observation)
    viterbi = zeros(N,T)
    first_obs = observation[1]
    backpointer = zeros(T)
    #initialization
    for i = 1:N
        viterbi[i,1] = pi[1,i]*e_matrix[first_obs,i]
    end
    #recursion
    for i = 2:T
        for k = 1:N
            previous = viterbi[i-1:i]
            max_ = max(previous[1]*t_matrix[1,k]*e_matrix[observation[i],k],previous[2]*t_matrix[2,k]*e_matrix[observation[i],k])
            location = argmax([previous[1]*t_matrix[1,k]*e_matrix[observation[i],k],previous[2]*t_matrix[2,k]*e_matrix[observation[i],k]])
            viterbi[k,i] = max_
            backpointer[i] = location
        end
    end
    #termination
    bestpathprob = maximum(viterbi[1:N,T])
    bestpathpointer = argmax(viterbi[1:N,T])
    bestpath = zeros(T)
    for i = 2:T
        bestpath[i-1] = backpointer[i]
    end
    bestpath[T] = bestpathpointer
    bestpath = join(bestpath,"->")
    return bestpath, bestpathprob
end
pi = [.8 .2]
e_matrix = [.2 .5; .4 .4; .4 .1]
t_matrix = [.6 .4;.5 .5]
observation = [3 1 3]
viterbi(pi,observation,t_matrix, e_matrix)
#("1.0->2.0->1.0", 0.007680000000000003)
