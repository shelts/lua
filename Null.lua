arg = { ... }

assert(#arg == 6, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed")
argSeed = 34086709
prng = DSFMT.create(argSeed)
-- npa = new parameter arrangement --

evolveTime       = tonumber(arg[1])
reverseOrbitTime = tonumber(arg[1]) / tonumber(arg[2])
rscale_l         = tonumber(arg[3])
rscale_d         = tonumber(arg[4])
mass_l           = tonumber(arg[5])
mass_d           = tonumber(arg[6])


evolveTime       = 3.87427734322731 
reverseOrbitTime = evolveTime / 1.01387196544634
rscale_l         = 1.29523584391164
rscale_d         = 2.99564367318775
mass_l           = 26.1017521350673
mass_d           = 131.258621913187
print(string.format("%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n",evolveTime,reverseOrbitTime,rscale_l, rscale_d,mass_l, mass_d))

evolveTime       = 0.000000001
reverseOrbitTime = evolveTime / 1.01387196544634
model1Bodies = 20000
totalBodies = model1Bodies

nbodyLikelihoodMethod = "EMD"
nbodyMinVersion = "1.62"


function makePotential()
   
return  nil
end

function get_timestep()
        --Mass of a single dark matter sphere enclosed within light rscale
    mass_enc_d = mass_d * (rscale_l)^3 * ( (sqr(rscale_l)+ sqr(rscale_d) ) )^(-3.0/2.0)

    --Mass of a single light matter sphere enclosed within dark rscale
    mass_enc_l = mass_l * (rscale_d)^3 * ( (sqr(rscale_l)+ sqr(rscale_d) ) )^(-3.0/2.0)

    s1 = cube(rscale_l) / (mass_enc_d + mass_l)
    s2 = cube(rscale_d) / (mass_enc_l + mass_d)
    
    --return the smaller time step
    if(s1 < s2) then
        s = s1
    else
        s = s2
    end
    
    -- I did it this way so there was only one place to change the time step. 
    t = (1/100) * sqrt( pi_4_3 * s)
--     print(string.format("%.15f \n", cube( 123456.78)))
--     print(string.format("%.15f \n", ( 123456.78)^(3)))
    return t
end


function makeContext()
   soften_length = (mass_l * rscale_l + mass_d * rscale_d) / (mass_d + mass_l)
   print(string.format("timestep = %.15f  eps2 = %.15f\n",get_timestep(), calculateEps2(totalBodies, soften_length)))
   return NBodyCtx.create{
      timeEvolve = evolveTime,
      timestep   = get_timestep(),
      eps2       = calculateEps2(totalBodies, soften_length),
      criterion  = "NewCriterion",
      useQuad    = true,
      theta      = 1.0
   }
end


-- Also required
function makeBodies(ctx, potential)
    local firstModel
    local finalPosition, finalVelocity = Vector.create(0, 0, 0), Vector.create(0, 0, 0)
    firstModel = predefinedModels.isotropic{
        nbody       = model1Bodies,
        prng        = prng,
        position    = finalPosition,
        velocity    = finalVelocity,
        mass1       = mass_l, 
        mass2       = mass_d,
        scaleRadius1 = rscale_l,
        scaleRadius2 = rscale_d,
        ignore      = true
    }
 
    
return firstModel
end

function makeHistogram()
    return HistogramParams.create{
     --Orphan Stream coordinate transformation angles
     phi = 128.79,
     theta = 54.39,
     psi = 90.70,
     lambdaStart = -75,
     lambdaEnd = 50,
     lambdaBins = 100,
     betaStart = -40,
     betaEnd = 40,
     betaBins = 1
}
end
