arg = { ... }

assert(#arg == 6, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed")

prng = DSFMT.create(argSeed)
-- npa = new parameter arrangement again --

evolveTime       = tonumber(arg[1])
reverseOrbitTime = tonumber(arg[1]) / tonumber(arg[2])
rscale_l         = tonumber(arg[3])
light_r_ratio    = tonumber(arg[4])
mass_l           = tonumber(arg[5])
light_mass_ratio = tonumber(arg[6])

model1Bodies = 20000
totalBodies = model1Bodies

nbodyLikelihoodMethod = "EMD"
nbodyMinVersion = "1.54"

mass_d   = mass_l / light_mass_ratio
rscale_d = rscale_l / light_r_ratio

-- print(evolveTime)
-- print(mass_d, mass_l)
function makePotential()
   return  Potential.create{
      spherical = Spherical.spherical{ mass  = 1.52954402e5, scale = 0.7 },
      disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
      halo      = Halo.logarithmic{ vhalo = 73, scaleLength = 12.0, flattenZ = 1.0 }
   }
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
    return t
end


function makeContext()
   soften_length = (mass_l*rscale_l + mass_d*rscale_d)/(mass_d+mass_l)
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
-- position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6))
function makeBodies(ctx, potential)
  local firstModel
  local finalPosition, finalVelocity = reverseOrbit{
      potential = potential,
      position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6)),
      velocity  = Vector.create(-156, 79, 107),
      tstop     = reverseOrbitTime,
      dt        = ctx.timestep / 10.0
  }
  
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
     lambdaStart = -150,
     lambdaEnd = 150,
     lambdaBins = 50,
     betaStart = -40,
     betaEnd = 40,
     betaBins = 1
}
end


