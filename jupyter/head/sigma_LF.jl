# at LF https://itis.swiss/virtual-population/tissue-properties/database/low-frequency-conductivity/
sigma=Dict(
	1=>3.48e-1, # white matter
 11=>3.48e-1, # white matter (tumor zone)
	2=>4.19e-1, # grey matter
	3=>1.88,    # CSF
	5=>1.48e-1, # skin
	6=>3.40e-1, # eye (mean value)
	7=>1.79e-2, # bone
	8=>1.8e-1,  # bone marrow (yellow)
	9=>6.62e-1, # blood
 10=>4.61e-1  # muscle
)
src_j=Dict(
	1=>0.0,	# white matter
    11=>0.0,	# white matter (tumor zone)
	2=>0.0,	# grey matter
	3=>2.0, # CSF
	5=>0.0, # skin
	6=>1.0, # eye (mean value)
	7=>0.0, # bone
	8=>0.0, # bone marrow (yellow)
	9=>0.0, # blood
 10=>0.0  # muscle
)

name = Dict( #
   1  => "Brain (White Matter)",
  11  => "Tumor",
   2  => "Brain (Grey Matter)",
   3  => "Cerebrospinal Fluid",
   5  => "Skin (Dry)",
   6  => "Eye (Cornea)",
   7  => "Bone (Cortical)",
   8  => "Bone Marrow (Yellow)",
   9  => "Blood",
  10  => "Muscles"
)
