arg = { ... }

assert(#arg == 6, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed")

prng = DSFMT.create(argSeed)
-- print('seed=', argSeed)
evolveTime       = tonumber(arg[1])
reverseOrbitTime = tonumber(arg[1]) / tonumber(arg[2])
rscale_l         = tonumber(arg[3])
light_r_ratio    = tonumber(arg[4])
mass_l           = tonumber(arg[5])
light_mass_ratio = tonumber(arg[6])
-- print(evolveTime)
model1Bodies = 1

totalBodies = model1Bodies

nbodyLikelihoodMethod = "EMD"
nbodyMinVersion = "1.50"

dwarfMass = mass_l / light_mass_ratio
rscale_t  = rscale_l / light_r_ratio
rscale_d  = rscale_t *  (1.0 - light_r_ratio)
mass_d    = dwarfMass * (1.0 - light_mass_ratio)

-- rscale_t  = rscale_l / light_r_ratio
-- rscale_d  = rscale_t *  (1.0 - light_r_ratio)
-- mass_d   = 0.0 --dwarfMass * (1.0 - light_mass_ratio)
-- print(mass_d, mass_l)


function makePotential()
   
return  nil
end

function get_timestep()
        --Mass of a single dark matter sphere enclosed within light rscale
    mass_enc_d = mass_d * (rscale_l)^3 * ( (sqr(rscale_l)+ sqr(rscale_d) ) )^(-3.0/2.0)

    --Mass of a single light matter sphere enclosed within dark rscale
    mass_enc_l = mass_l * (rscale_d)^3* ( (sqr(rscale_l)+ sqr(rscale_d) ) )^(-3.0/2.0)

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
    return t
end


function makeContext()
   soften_length = (mass_l * rscale_l + mass_d * rscale_d) / (mass_d + mass_l)
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
     phi = 128.79,
     theta = 54.39,
     psi = 90.70,
     lambdaStart = 0.0,
     lambdaEnd = 300,
     lambdaBins = 200,
     betaStart = 0.0,
     betaEnd = -200,
     betaBins = 1
}
end
