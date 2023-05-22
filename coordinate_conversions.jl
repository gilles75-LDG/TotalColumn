
lat = 42.2
lon = -79.4
θ₀ = 90 - lat
ϕ₀ = lon

function radial_from_perpPlane(x,y,r=1,θ₀=0,ϕₒ=0)
    #returns the change in relative radial coordinates for an x,y position on a perpendicular tangential plane (N/S = y, E/W = x) 
    #returns r,θ,ϕ (r, 90-lat, lon)
    return (r*secd(sqrt(x^2 + y^2)), θ₀+atand(y/r), ϕₒ+atand(x/r) )
end

radial_from_perpPlane(100,100,6000)

#TODO planar projection from ASZA and AZIM

