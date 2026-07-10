function barea( b )
    #println( maximum(b[1])-minimum(b[1]), ' ', maximum(b[2])-minimum(b[2]) )
    return (maximum(b[1])-minimum(b[1]))*(maximum(b[2])-minimum(b[2]))
end

function bbintersect( b1, b2 )
    x1l,x1h,y1l,y1h = minimum(b1[1]), maximum(b1[1]), minimum(b1[2]), maximum(b1[2])
    x2l,x2h,y2l,y2h = minimum(b2[1]), maximum(b2[1]), minimum(b2[2]), maximum(b2[2])
    xover = maximum([0,minimum([x1h,x2h])-maximum([x1l,x2l])])
    yover = maximum([0,minimum([y1h,y2h])-maximum([y1l,y2l])])
    return xover,yover
end
    
function iandu( b1, b2 )
    xover,yover = bbintersect( b1, b2 )
    inte = xover*yover
    xunion = maximum([b1[1];b2[1]])-minimum([b1[1];b2[1]])
    yunion = maximum([b1[2];b2[2]])-minimum([b1[2];b2[2]])
    unio = (maximum(b1[1])-minimum(b1[1]))*(maximum(b1[2])-minimum(b1[2])) + (maximum(b2[1])-minimum(b2[1]))*(maximum(b2[2])-minimum(b2[2])) - inte
    return inte,unio,xover/xunion,yover/yunion
end

function bunion( b1, b2 )
    xl,xh,yl,yh = minimum([b1[1];b2[1]]), maximum([b1[1];b2[1]]), minimum([b1[2];b2[2]]), maximum([b1[2];b2[2]])
    return [[xl,xh,xh,xl,xl],[yl,yl,yh,yh,yl]]
end

function compactbblist( bb; iou=0.45, xyiou=(0.05,0.9), maxit=10 )
    aa= []
    append!(aa, bb)
    compacted= true
    it= 0
    while compacted && it < maxit
        compacted = false
        for i=1:length(aa)
            if aa[i][1][1] >= 0
                for j=i+1:length(aa)
                    if aa[j][1][1] >= 0
                        is,u,xi,yi=iandu(aa[i],aa[j])
                        ai=barea(aa[i])
                        aj=barea(aa[j])
                        mina=min(ai,aj)
                        if ( is/mina > iou ) || (is/mina > xyiou[1] && (xi > xyiou[2] || yi > xyiou[2] || xi > xyiou[1] || yi > xyiou[1] ))
                            cc = bunion( aa[i], aa[j])
                            aa[i] = cc
                            aa[j] = [[-2,-1,-1,-2],[-2,-2,-1,-1]]
                            compacted = true
                        end
                    end
                end
            end
        end
        it += 1
    end

    cc=[]
    for b in aa
        if b[1][1] >= 0
            push!(cc,b)
        end
    end
    return cc
end

function saveBB( bb, dx, dy, f )
    open( f, "w" ) do file
        for b in bb
            xl,xh,yl,yh=minimum(b[1]),maximum(b[1]),minimum(b[2]),maximum(b[2])
            xb,yb,wb,hb = 0.5*(xl+xh)/dx,0.5*(yl+yh)/dy,(xh-xl)/dx,(yh-yl)/dy
            write( file, "0 $xb $yb $wb $hb\n")
        end
    end
end

function saveFLD( x2, y2, z2, E2, f )
    open(f, "w") do file
        for i=1:size(x2)[1]
            vx = x2[i]
            for j = 1:size(y2)[1]
                vy = y2[j]
                for k = 1:size(z2)[1]
                    vz = z2[k]
                    vE = E2[i,j,k]
                    write( file, "$vx $vy $vz $vE\n")
                end
            end
        end
    end
end

function pkl2txt( f, thresh, maxit )
    mesh,V,err,bp = uziom_adaptive(f,maxit)
    x=mesh["x"]
    nop=mesh["v"]
    bb,E,mE=hotspots(x,nop,V,thresh,gradVIn3D4NElem)
    
    dx,dy=maximum(x[:,1])-minimum(x[:,1]),maximum(x[:,2])-minimum(x[:,2])
    bbb=compactbblist(bb, maxit=3)
    
    saveBB( bbb, dx, dy, replace(f,"pkl" => "txt") )
    
    #x2,y2,z2,E2,i2,errv=fastScan3DCnst(x, nop, mE, 0)

    #saveFLD( x2, y2, z2, E2, replace(f,"pkl" => "fld") )

    return Dict("mesh"=>mesh,"V"=>V,"err"=>err,"BPData"=>bp)
end
