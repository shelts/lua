
arg = { ... }
-- /* Copyright (c) 2016 Siddhartha Shelton */
assert(#arg == 4, "Expected 4 arguments")
assert(argSeed ~= nil, "Expected seed")
argSeed = 34086709
prng = DSFMT.create(argSeed)

evolveTime       = tonumber(arg[1])
reverseOrbitTime = evolveTime / tonumber(arg[2])

print(evolveTime)

r0  = tonumber(arg[3])
r0 = 0.2
dwarfMass  = tonumber(arg[4])
dwarfMass = 60
model1Bodies = 20000
totalBodies = model1Bodies

function makePotential()
   return  Potential.create{
      spherical = Spherical.spherical{ mass  = 1.52954402e5, scale = 0.7 },
      disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
      halo      = Halo.logarithmic{ vhalo = 73, scaleLength = 12.0, flattenZ = 1.0 }
   }
end



encMass = plummerTimestepIntegral(r0, dwarfMass, 1e-7)

-- This is also required
function makeContext()
   return NBodyCtx.create{
      timeEvolve = evolveTime,
      timestep   = sqr(1/10.0) * sqrt((pi_4_3 * cube(r0)) / (dwarfMass)),
      eps2       = calculateEps2(totalBodies, r0),
      criterion  =  "NewCriterion",
      useQuad    = true,
      theta      = 1.0
   }
end

-- Also required
function makeBodies(ctx, potential)
   local firstModel
  local finalPosition, finalVelocity = reverseOrbit{
      potential = potential,
      position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6)),
      velocity  = Vector.create(-156, 79, 107),
      tstop     = reverseOrbitTime,
      dt        = ctx.timestep / 10.0
  }

   firstModel = predefinedModels.plummer{
      nbody       = model1Bodies,
      prng        = prng,
      position    = finalPosition,
      velocity    = finalVelocity,
      mass        = dwarfMass,
      scaleRadius = r0,
      ignore      = false
   }

   return firstModel
end

function makeHistogram()
    return HistogramParams.create{
     phi = 128.79,
     theta = 54.39,
     psi = 90.70,
     lambdaStart = -180,
     lambdaEnd = 180,
     lambdaBins = 50,
     betaStart = -180,
     betaEnd = 180,
     betaBins = 1
}
end