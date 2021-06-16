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

#returns stationary distribution given transition matrix through eigendecomposition
#using QuantEcon
function s_stationary(t_matrix)
    mc = MarkovChain(t_matrix)
    answer = stationary_distributions(mc)
    #answer = split(answer)
    return answer
end
B = [.2 .6 .2 ; .3 0 .7 ; .5 0 .5 ]
s_stationary(B)  

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
