# q is the order of hmm
#returns most likely sequence from observed values
function viterbi(pi, obs, tm, em,q)
    N = length(pi)
    T = length(obs)
    viterbi = zeros(N,T)
    first_obs = obs[1]
    backpointer = zeros(T)
    #initialization
    
    for i = 1:N
        viterbi[i,1] = pi[1,i]*em[first_obs][i]
    end
    
    for i = 2:q	
    	for k = 1:N
    		prev = viterbi[1:N, i-1]
            max_ = max(prev[1]*tm[1,k]*em[obs[i]][k],prev[2]*tm[2,k]*em[obs[i]][k])
            location = argmax([prev[1]*tm[1,k]*em[obs[i]][k],prev[2]*tm[2,k]*em[obs[i]][k]])
            viterbi[k,i] = max_
            backpointer[i] = location
        end
    end
    		
    #      ????????
    #calculations
    for i = q+1:T
        for k = 1:N
            prev = viterbi[1:N, i-1]
            #println(prev)
            max_ = max(prev[1]*tm[1,k]*em[obs[i]][k],prev[2]*tm[2,k]*em[obs[i]][k])
            location = argmax([prev[1]*tm[1,k]*em[obs[i]][k],prev[2]*tm[2,k]*em[obs[i]][k]])
            #location = argmax(max_)
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
    #bestpath = join(bestpath,"->")
    return bestpath 

end

pi = [1 0]
em = Dict([("A",[.1 .4]),
			("C",[.4 .1]),
			("G",[.4 .1]),
			("T",[.1 .4])])
#[exon intron]
tm = [.99 .01;.02 .98]
#=
EE	EI
IE	II
=#
obs = split("GGCCCCCGCGCGCGCGGGGGGCCCATATATATATATTTTTTAAGCGCGGCGCGCGGCGCGCATATATATATATA","")	
#			 EEEEEEEEEEEEEEEEEEEEEEEEIIIIIIIIIIIIIIIIIIIIIIIIIEEEEEEEEEEEEIIIIIIIIIIIII	o
for element in viterbi(pi,obs,tm, em,3)
	if element == 2.0 print('I')
	elseif element == 1.0 print('E') end
end
println()