use Module::Build;
# See perldoc Module::Build for details of how this works

Module::Build->new
    ( module_name     => 'PerlIO::via::SeqIO',
      license         => 'perl',
      perl            => 5.8,
      requires        => {
                          'Bio::Root::Version'=> '1.006',
			  'IO::String'     => 0,
			  'IO::Seekable'   => 0
                         },
      build_requires  => {
	                  'Module::Build'  => 0,
                          'Test::More'     => 0,
	                  'File::Temp'     => 0,
			  'File::Spec'     => 0,
                         },
      create_makefile => 'passthrough'
    )->create_build_script;
