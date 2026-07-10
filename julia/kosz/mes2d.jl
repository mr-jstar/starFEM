include( "mes.jl" )

x_min,y_min = 0,0
x_max = 1.0
y_max = 0.7
nx = 50
ny = 20

x,nop = box2D( x_min, y_min, x_max, y_max, nx, ny )

f=plotMesh( x, nop )

print(size(x), size(nop))

@show f
