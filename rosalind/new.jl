# trying to make own fasta thing

function fasta_decoder(file)
	file = readlines(open(file, "r"))
	id = []
	seqs = []
	i = 1
	while i <= length(file)
		seq = ""
		if file[i][1] == '>' 
			push!(id, file[i])
		end
		i += 1
		while file[i][1] != '>'
			seq = seq*file[i]
			i += 1
			if i > length(file) break end
		end
		push!(seqs, seq)
	end
	return(id, seqs)
end
println(fasta_decoder("test.fa"))