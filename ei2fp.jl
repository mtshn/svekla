using SparseArrays
using Statistics
using MLUtils
using ProgressBars
using Flux
using Flux.Data: DataLoader
using Random
using DelimitedFiles
using CUDA



LENGTH_FP = 1024
LENGTH_SPECTRUM = 800

function create_sparse_vectors_from_java_string(s,l)
    strs = split(s)
    n_nonzeros = Int(round(length(strs)/2))

    x1 = 1
    result = spzeros(Float32, Int32, l)

    for i = 1:n_nonzeros
        x = parse(Int32, strs[x1]) + 1
        x1 = x1 + 1
        y = parse(Float32, strs[x1]) 
        x1 = x1 + 1
        result[x] = y
    end

    return result
end

function load_tabular_dataset(filename1,filename2,l1,l2,scale1,scale2)
    lines1 = readlines(filename1)
    lines_no_empty1 = String[]
    for l in lines1
        if strip(l) != ""
            push!(lines_no_empty1, l)
        end
    end
    lines2 = readlines(filename2)
    lines_no_empty2 = String[]
    for l in lines2
        if strip(l) != ""
            push!(lines_no_empty2, l)
        end
    end
    result = Array{Tuple{SparseVector{Float32},SparseVector{Float32}},1}(undef, length(lines_no_empty1))
    for i in eachindex(lines_no_empty1)
        line1 = lines_no_empty1[i]
        line2 = lines_no_empty2[i]
        g1 = create_sparse_vectors_from_java_string(line1,l1)
        g2 = create_sparse_vectors_from_java_string(line2,l2)
        result[i] = (g1*scale1,g2*scale2)
    end
    return result
end

function x_y_batch(sparse_data_batch)
    t = sparse_data_batch    
    n = length(sparse_data_batch)
    result1 = zeros(Float32, length(t[1][1]), n)
    result2 = zeros(Float32, length(t[1][2]), n)
    for i = 1:n
        for j = 1:length(t[1][1])
            result1[j, i] = t[i][1][j]
        end
        for j = 1:length(t[1][2])
            result2[j, i] = t[i][2][j]
        end
    end
    return result1, result2
end


function loss_cross_entropy(x, y, mdl)
    y1 = mdl(x)
    l = Flux.logitbinarycrossentropy(y1, y)
    return l
end

function train(i, train_set,val_set,parameters)
    for g in ProgressBar(train_set)
        x, y = x_y_batch(g)
        x = x |> Flux.gpu
        y = y |> Flux.gpu
        gs = Flux.gradient(parameters) do
            loss_cross_entropy(x, y, mdl)
        end
        Flux.Optimise.update!(opt, parameters, gs)
    end
    #validate(train_set)
    return validate(val_set,mdl)
end

function predict(dataloader, mdl)
    y1 = Array{Float32,2}[]
    for g in dataloader
        x, y = x_y_batch(g)
        x = x |> Flux.gpu
        y = mdl(x) |> Flux.cpu
        push!(y1, y)
    end
    return (reduce(hcat, y1))
end


function labels(dataloader)
    y1 = Array{Float32,2}[]
    for g in dataloader
        _, y = x_y_batch(g)
        push!(y1, y)
    end
    return (reduce(hcat, y1))
end


function validate(dataloader,  mdl)
    y2 = predict(dataloader,mdl)
    y1 = labels(dataloader)

    l = Flux.logitbinarycrossentropy(y2, y1)

    return l
end

struct EI2FP
    dense1
    dense2
    dropout1
    dense3
    dropout2
    dense4
    batchnorm
    dense5
end

silu(x)=x*sigmoid(x)

Flux.@functor EI2FP

function EI2FP()
    EI2FP(
        Dense(LENGTH_SPECTRUM  => 4096, silu),
        Dense(4096 => 4096, silu),
        Dropout(0.5),
        Dense(4096 => 2048, silu),
        Dropout(0.5),
        Dense(2048 => 2048, silu),
        BatchNorm(2048),
        Dense(2048 => LENGTH_FP, identity),
    )
end

function (EI2FP::EI2FP)(x)
    x = EI2FP.dropout2(EI2FP.dense3(EI2FP.dropout1(EI2FP.dense2(EI2FP.dense1(x)))))
    x = EI2FP.batchnorm(EI2FP.dense4(x))
    x = EI2FP.dense5(x)
    return x
end

##########


batch_size = 256

train_val_array = load_tabular_dataset("../train_ei2fp/spectra.txt","../train_ei2fp/fp.txt",LENGTH_SPECTRUM,LENGTH_FP,0.001,1)
sample = randsubseq(1:size(train_val_array, 1), 0.1)
not_sample = [i for i in 1:size(train_val_array, 1) if isempty(searchsorted(sample, i))]

train_array = train_val_array[not_sample]
val_array = train_val_array[sample]

train_set = DataLoader(shuffle(train_array), batchsize=batch_size)
val_set = DataLoader(val_array, batchsize=batch_size)


mdl = EI2FP() |> Flux.gpu
ps = Flux.params(mdl) |> Flux.gpu
opt = ADAM(0.001)
best_loss = 1000000
best_params = collect(Flux.params(cpu(mdl)))
for i = 1:175
    if i%50 ==0
        opt.eta=opt.eta/2
    end
    loss =  train(i, train_set, val_set,ps)
    println(i)
    println(loss)
    if loss < best_loss
        best_loss = loss
        best_params = deepcopy(collect(Flux.params(cpu(mdl))))
    end
end
Flux.loadparams!(mdl, best_params)

mdl = mdl |> Flux.cpu
folder = "../train_ei2fp/params/"
writedlm(string(folder,"dense1_weight.txt"),mdl.dense1.weight)
writedlm(string(folder,"dense2_weight.txt"),mdl.dense2.weight)
writedlm(string(folder,"dense3_weight.txt"),mdl.dense3.weight)
writedlm(string(folder,"dense4_weight.txt"),mdl.dense4.weight)
writedlm(string(folder,"dense5_weight.txt"),mdl.dense5.weight)
writedlm(string(folder,"dense1_bias.txt"),mdl.dense1.bias)
writedlm(string(folder,"dense2_bias.txt"),mdl.dense2.bias)
writedlm(string(folder,"dense3_bias.txt"),mdl.dense3.bias)
writedlm(string(folder,"dense4_bias.txt"),mdl.dense4.bias)
writedlm(string(folder,"dense5_bias.txt"),mdl.dense5.bias)
writedlm(string(folder,"batchnorm_mu.txt"),mdl.batchnorm.μ)
writedlm(string(folder,"batchnorm_beta.txt"),mdl.batchnorm.β)
writedlm(string(folder,"batchnorm_gamma.txt"),mdl.batchnorm.γ)
writedlm(string(folder,"batchnorm_sigma2.txt"),mdl.batchnorm.σ²)
