arg = { ... }

assert(#arg == 6, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed")
argSeed = 34086709
prng = DSFMT.create(argSeed)

evolveTime       = tonumber(arg[1])
rev_ratio        = tonumber(arg[2])
rscale_l         = tonumber(arg[3])
light_r_ratio    = tonumber(arg[4])
mass_l           = tonumber(arg[5])
light_mass_ratio = tonumber(arg[6])

function round(num, places)
  local mult = 10.0^(places)
  return floor(num * mult + 0.5) / mult
end

dec = 9.0
evolveTime       = round( evolveTime,       dec )
rev_ratio        = round( rev_ratio,        dec )
rscale_l         = round( rscale_l,         dec )
light_r_ratio    = round( light_r_ratio,    dec )
mass_l           = round( mass_l,           dec )
light_mass_ratio = round( light_mass_ratio, dec )
model1Bodies = 20000
totalBodies = model1Bodies

nbodyLikelihoodMethod = "EMD"
nbodyMinVersion = "1.62"

revOrbTime = evolveTime / rev_ratio
dwarfMass = mass_l / light_mass_ratio
rscale_t  = rscale_l / light_r_ratio
rscale_d  = rscale_t *  (1.0 - light_r_ratio)
mass_d    = dwarfMass * (1.0 - light_mass_ratio)

charles_run = true
print_reverse_orbit = true

if(charles_run == true) then
    evolveTime = 2.0
    revOrbTime = 2.0
    
    dwarfMass = 5e6 / 222288.47
    mass_ratio = 1/10
    
    mass_d = dwarfMass / (mass_ratio + 1)
    mass_l = dwarfMass - mass_d
    mass_d = 2.0 * mass_d
    rscale_d = .25
    r_ratio = 1.0/5.0
    rscale_l = rscale_d * r_ratio
    
    mass_d = 5e6 / 222288.47
    rscale_d = 0.125
    mass_l = 1e5 / 222288.47
    rscale_l = 0.01

    print('forward time=', evolveTime, '\nrev time=',  revOrbTime)
    print('mass_l sim=', mass_l, '\nmass_d sim=', mass_d)
    print('light mass solar=', mass_l * 222288.47, '\ndark mass solar=', mass_d * 222288.47)
    print('total mass solar= ', (mass_d + mass_l) * 222288.47)
    print('rl = ', rscale_l, 'rd = ', rscale_d)
end

function makePotential()
   return  Potential.create{
      spherical = Spherical.spherical{ mass  = 1.52954402e5, scale = 0.7 },
      disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
      halo      = Halo.nfw{ vhalo = 73, scaleLength = 22.250}
   }
end


function get_timestep()
    --Mass of a single dark matter sphere enclosed within light rscale
    mass_enc_d = mass_d * (rscale_l)^3 * ( (rscale_l)^2 + (rscale_d)^2  )^(-3.0/2.0)

    --Mass of a single light matter sphere enclosed within dark rscale
    mass_enc_l = mass_l * (rscale_d)^3 * ( (rscale_l)^2 + (rscale_d)^2  )^(-3.0/2.0)

    s1 = (rscale_l)^3 / (mass_enc_d + mass_l)
    s2 = (rscale_d)^3 / (mass_enc_l + mass_d)
    
    --return the smaller time step
    if(s1 < s2) then
        s = s1
    else
        s = s2
    end
    
    -- I did it this way so there was only one place to change the time step. 
    t = (1 / 100.0) * ( pi_4_3 * s)^(1.0/2.0)
    
--     t = sqr(1/10.0) * sqrt((pi_4_3 * cube(rscale_d)) / (dwarfMass))
    print('timestep ', t)
    
    return t
end


function makeContext()
   soften_length = (mass_l * rscale_l + mass_d  * rscale_d) / (mass_d + mass_l)
   return NBodyCtx.create{
      timeEvolve = evolveTime,
      timestep   = get_timestep(),
      eps2       = calculateEps2(totalBodies, soften_length ),
      criterion  = "NewCriterion",
      useQuad    = true,
      theta      = 1.0
   }
end

if(charles_run == true) then
    soften_length = (mass_l * rscale_l + mass_d  * rscale_d) / (mass_d + mass_l)
    print('soften_length ', calculateEps2(totalBodies, soften_length ))
end
-- Also required
-- for orphan: lbr in sun centered
-- position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6))
-- velocity  = Vector.create(-156, 79, 107),
function makeBodies(ctx, potential)
  local firstModel
  local finalPosition, finalVelocity = reverseOrbit{
      potential = potential,
      position  = lbrToCartesian(ctx, Vector.create(45, 46.93, 11.87)),
      velocity  = Vector.create(-122.78, 157.32, 64.90),
      tstop     = revOrbTime,
      dt        = ctx.timestep / 10.0
    }
  
    
  if(print_reverse_orbit == true) then
    local placeholderPos, placeholderVel = PrintReverseOrbit{
        potential = potential,
        position  = lbrToCartesian(ctx, Vector.create(45, 46.93, 11.87)),
        velocity  = Vector.create(-122.78, 157.32, 64.90),
        tstop     = .14,
        tstopf    = .20,
        dt        = ctx.timestep / 10.0
    }
    print('Printing reverse orbit')
  end
  
  if(charles_run == true) then
    print(lbrToCartesian(ctx, Vector.create(45, 46.93, 11.87)), Vector.create(-122.78, 157.32, 64.90))
    print(finalPosition, finalVelocity)
  end
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
--     for i, v in ipairs(firstModel)
--     do
--         print(firstModel[i].position.x, firstModel[i].position.y, firstModel[i].position.z, firstModel[i].velocity.x, firstModel[i].velocity.y, firstModel[i].velocity.z)
--         print(firstModel[i].mass)
--     end
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
     lambdaBins = 100,
     betaStart = -100,
     betaEnd = 100,
     betaBins = 1
}
end


