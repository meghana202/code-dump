function perm(n,m)
	list = []
	nums = 1:n
	if m == 1
		for i = 1:n push!(list, string(i)) end
		return list
	elseif m > 1
		for i = 1:length(perm(n,m-1)) 
			for j in nums
				if !occursin(string(j), perm(n, m-1)[i])
					push!(list, perm(n,m-1)[i]*string(j)) 
				end 
			end
		end 
		return list
	end			
end
println(perm(5,5))


#=
function output_formatting(n)
	println(factorial(n))
	for i = 1:length(perm(n,n)) println(perm(n,n)[i]) end
end
output_formatting(5)
=#