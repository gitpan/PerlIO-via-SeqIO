use Module::Build;
# See perldoc Module::Build for details of how this works

Module::Build->new
    ( module_name     => 'PerlIO::via::SeqIO',
      license         => 'perl',
      requires        => {
	                  perl            => 5.010,
                          'Bio::Root::Version'=> '1.006',
			  'IO::String'     => 0,
			  'IO::Seekable'   => 0,
			  'PerlIO::Util'   => 0,
			  'Scalar::Util'   => 1.19,
			  'Symbol'         => 0
                         },
      recommends      => {
	                  'PerlIO::via::gzip' => 0.031
                         },
      build_requires  => {
	                  'Module::Build'  => 0,
                          'Test::More'     => 0,
	                  'File::Temp'     => 0,
			  'File::Spec'     => 0,
                         },
      create_makefile => 'passthrough'
    )->create_build_script;
