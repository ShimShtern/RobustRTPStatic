using Plots
β=[1e-7,5e-7,1e-6,1e-5,2e-5,3e-5,4e-5,5e-5,6e-5,7e-5,8e-5,9e-5,1e-4]
g = [78.50927224871357,78.50836787666066,
78.50597598977697,78.44034671247854,
78.38903191749267,78.33675456561923,78.29939210564523,
78.2657747775239,78.23083117779953,78.19117772317935,78.15246877395913,78.10969241038319,78.06828956968708]
nonZeroDev = [49610,46338,43238,23959,19020,16137,14477,
13306,12344,11512,10772,10070, 9458]

DVRHS = 13601
#l = @layout [a ; b ]
#gr()
#p1=plot(β,g,seriestype = :scatter,xlabel="β",ylabel="g",label="")

#p2=plot(β,nonZeroDev,seriestype = :scatter,xlabel="β",ylabel="||δ||_0",label="")
p1=plot(β,g,xlabel="β",ylabel="g")
p2=plot(β,nonZeroDev,xlabel="β",ylabel="||δ||_0",label="DV LHS")
plot!([0,maximum(β)],[DVRHS,DVRHS],label="DV RHS")

plot(p1,p2,layout=2)
