# Configuration file for peptidesim.

#------------------------------------------------------------------------------
# Configurable configuration
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# LoggingConfigurable configuration
#------------------------------------------------------------------------------

# A parent class for Configurables that log.
# 
# Subclasses have a log trait, and the default behavior is to get the logger
# from the currently running Application.

#------------------------------------------------------------------------------
# SingletonConfigurable configuration
#------------------------------------------------------------------------------

# A configurable that only allows one instance.
# 
# This class is for classes that should only have one instance of itself or
# *any* subclass. To create and retrieve such a class use the
# :meth:`SingletonConfigurable.instance` method.

#------------------------------------------------------------------------------
# Application configuration
#------------------------------------------------------------------------------

# This is an application.

# The date format used by logging formatters for %(asctime)s
# c.Application.log_datefmt = '%Y-%m-%d %H:%M:%S'

# The Logging format template
# c.Application.log_format = '[%(name)s]%(highlevel)s %(message)s'

# Set the log level by value or name.
# c.Application.log_level = 30

#------------------------------------------------------------------------------
# PeptideSimConfigurator configuration
#------------------------------------------------------------------------------

# The config file to load
# c.PeptideSimConfigurator.config_file = 'peptidesim_config.py'

# Generate default config file
# c.PeptideSimConfigurator.generate_config = True

#------------------------------------------------------------------------------
# PeptideSim configuration
#------------------------------------------------------------------------------

# PeptideSim
# 
# Example -------
# 
# Here's an example for creating a ``PeptideSim`` object with the peptide AEAE
# using the default configuration, saving in the current directory. ::
# 
#     p = PeptideSim( dir_name = ".", seqs = ['AEAE'], counts = [1] )
# 
# Here's an example showing **one** AEAE peptide and **two** LGLG peptides,
# saving in the current directory. ::
# 
#     p = PeptideSim( dir_name = ".", seqs = ['AEAE', 'LGLG'], counts = [1,2])
# #counts in order of the list of peptides

# The config file to load
# c.PeptideSim.config_file = 'peptidesim_config.py'

# The name of the combined peptides without water
# c.PeptideSim.dry_packed_file = u'dry_mix.pdb'

# Files to carry over for each simulation
# c.PeptideSim.file_list = []

# The Gromcas structure file
# c.PeptideSim.gro = u'protein.gro'

# Initial PDB files for peptides
# c.PeptideSim.pdbfiles = []

# The density of the peptides in milligrams / milliliter
# c.PeptideSim.peptide_density = 0.2

# Barostat pressure. Ignored if not doing NPT
# c.PeptideSim.pressure = 0

# The name for the type of simulation job (e.g., NVE-equil-NVT-prod)
# c.PeptideSim.sim_name = u'peptidesim'

# Gromacs topology file
# c.PeptideSim.topol = u'topology.top'
