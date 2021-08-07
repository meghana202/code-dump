using ArgParse

cli = ArgParseSettings()
@add_arg_table cli begin
	"fasta"
		help = "fasta file, required"
		required = true
end
arg = parse_args(ARGS, cli)

#= if it was a first order thing
function uno(fasta,n)
	nt = Dict([('A', 0.0),('C', 0.0),('G', 0.0),('T', 0.0)])
	count_ = 0
	for line in fasta
		if line[1] != '>'
			#println(length(line))
			count_ += length(line)
			nt['A'] += (count(i->(i=='A'), line))
			nt['C'] += (count(i->(i=='C'), line))
			nt['G'] += (count(i->(i=='G'), line))
			nt['T'] += (count(i->(i=='T'), line))
		end
	end
	for (key,value) in nt
		nt[key] = value/count_ 
	end
	println(nt)
end

fasta = readlines(open(arg["fasta"], "r"))
uno(fasta,1)
=#


#ok now second order
function two(fasta, n)
	nt = Dict('A'=>0, 'C'=>0, 'G'=>0, 'T'=>0)
	println(nt)

end
fasta = readlines(open(arg["fasta"], "r"))
two(fasta, 2)