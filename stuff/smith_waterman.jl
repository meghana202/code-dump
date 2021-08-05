#match = 1, mismatch = -1, gap = -5
function sw()
	seq1 = ARGS[1]
	seq2 = ARGS[2]
	matrix_ = fill(0,length(seq2)+1,length(seq1)+1)
	trace_ = Matrix{Union{Char, Nothing}}(nothing,length(seq2)+1,length(seq1)+1)

	#first row and column
	for i = 1:length(seq1)
		trace_[1, i+1] = 'L'
	end
	for i = 1:length(seq2)
		trace_[i+1, 1] = 'U'
	end
	trace_[1,1] = '-' 
	
	#fill
	for i = 2:length(seq2)+1
		for j = 2:length(seq1)+1
			#x stores whether it was a match or mismatch/indel
			if seq1[j-1]==seq2[i-1] x = 1 else x = -1 end
			
			#filling matrix
			top = matrix_[i-1,j] - 5
			left = matrix_[i, j-1] - 5
			di = matrix_[i-1,j-1] + x
			if max(top, left, di) <= 0 matrix_[i,j] = 0
			else matrix_[i,j] = max(top, left, di) end
			
			#filling trace
			if max(top, left, di) == top trace_[i,j] = 'U'
			elseif max(top, left, di) == left trace_[i,j] = 'L'
			elseif max(top, left, di) == di trace_[i,j] = 'D' end
			
		end
	end

	#trace
	seqs = []
	max_ = maximum(matrix_)
	length_ = length(findall(x->x == max_, matrix_))
	for i = 1:length_
		cell_ = findall(x->x == max_, matrix_)[i]
		firstseq = ""
		secondseq = ""
		x = cell_[1]
		y = cell_[2]
		while true
			x = cell_[1]
			y = cell_[2]
			cell_ = [x, y]
			if matrix_[x,y] == 0 break end
			if trace_[x,y] == 'U'
				secondseq *= seq2[x-1]
				firstseq *= "-"
				cell_ = [x-1,y]
			elseif trace_[x,y] == 'L'
				firstseq *= seq1[y-1]
				secondseq *= "-"
				cell_ = [x, y-1]
			elseif trace_[x,y] == 'D'
				secondseq *= seq2[x-1]
				firstseq *= seq1[y-1]
				cell_ = [x-1, y-1]
			end
		end
		
		push!(seqs,(firstseq, secondseq))
	end
	
	#output/formatting
	longest = 0 
	for i = 1:length(seqs)
		if length(seqs[i][1]) > longest
			longest = length(seqs[i][1])
			global answer = seqs[i]
		end
	end
	for element in answer println(reverse(element)) end
end
sw()
