AUTOMAKE_OPTIONS = subdir-objects

#append flags
#ACLOCAL_AMFLAGS = -I m4 --install ${ACLOCAL_FLAGS}

lib_LIBRARIES = libmct.a
libmct_a_SOURCES =  mct/model.cpp mct/model.hpp \
				mct/model_options_parser.cpp mct/model_options_parser.hpp \
				mct/arg_reader.cpp mct/arg_reader.hpp \
				mct/arg_parser.cpp mct/arg_parser.hpp \
				mct/nodes_and_trees.cpp mct/nodes_and_trees.hpp \
				mct/model_set_config_builder.cpp mct/model_set_config_builder.hpp \
				mct/modelfactory.cpp mct/modelfactory.hpp \
				mct/modelset.cpp mct/modelset.hpp \
				mct/multi_loci_modelset_set.cpp mct/multi_loci_modelset_set.hpp \
				mct/mct_arg_factory.cpp mct/mct_arg_factory.hpp \
				mct/partition.cpp mct/partition.hpp \
				mct/partitioned_arg_sequence_generator.cpp \
				mct/partitioned_arg_sequence_generator.hpp \
				mct/mct_sample_sequence.cpp mct/mct_sample_sequence.hpp \
				mct/ind_site_sequence.cpp mct/ind_site_sequence.hpp \
				mct/model_data.hpp \
				mct/site_seq.hpp \
				mct/mct_model_data_types.hpp \
				mct/mct_locus_map.cpp mct/mct_locus_map.hpp \
				mct/genetic_pop_structure_analysable.cpp mct/genetic_pop_structure_analysable.hpp \
				mct/poly_data_sequence_maker.cpp mct/poly_data_sequence_maker.hpp \
				mct/poly_data_sequence_reader.cpp mct/poly_data_sequence_reader.hpp \
				mct/multi_loci_polytable.cpp mct/multi_loci_polytable.hpp \
				mct/multi_loci_polytable_set.cpp mct/multi_loci_polytable_set.hpp \
				mct/multi_loci_simdata_sampler.cpp mct/multi_loci_simdata_sampler.hpp \
				mct/multi_loci_polysites_sampler.cpp mct/multi_loci_polysites_sampler.hpp \
				mct/multi_loci_polytable_sampler.cpp mct/multi_loci_polytable_sampler.hpp \
				mct/mct_microsat_locus.cpp mct/mct_microsat_locus.hpp \
				mct/multi_loci_microsat.cpp mct/multi_loci_microsat.hpp \
				mct/multi_loci_microsat_set.cpp mct/multi_loci_microsat_set.hpp \
				mct/multi_loci_microsat_sampler.cpp mct/multi_loci_microsat_sampler.hpp \
				mct/population_labeled_arg.cpp mct/population_labeled_arg.hpp \
				mct/libseq.cpp mct/libseq.hpp \
				mct/arg_ms_converter.cpp mct/arg_ms_converter.hpp \
				mct/prng.hpp \
				mct/prng_gsl.cpp mct/prng_gsl.hpp \
				mct/summary_statistic.cpp mct/summary_statistic.hpp \
				mct/summary_statistic_set.cpp mct/summary_statistic_set.hpp \
				mct/descriptive_stats.cpp mct/descriptive_stats.hpp \
				mct/ahabc_values.cpp mct/ahabc_values.hpp \
				mct/utilities.cpp mct/utilities.hpp \
				aabc/rmv.cpp aabc/rmv.h \
				aabc/data_object.hpp \
				aabc/data_object_set.cpp aabc/data_object_set.hpp \
				aabc/generic_data_object_set.cpp aabc/generic_data_object_set.hpp \
				aabc/generic_data_object_sampler.cpp aabc/generic_data_object_sampler.hpp \
				aabc/parameter_particle_validator_generic_data_object.cpp aabc/\
				aabc/parameter_particle_validator_generic_data_object.hpp \
				aabc/parameter_particle.cpp aabc/parameter_particle.hpp \
				aabc/parameter_particle_set.cpp aabc/parameter_particle_set.hpp \
				aabc/parameter_particle_set_analysis.cpp aabc/parameter_particle_set_analysis.hpp \
				aabc/parameter_particle_validation_result.cpp aabc/parameter_particle_validation_result.hpp \
				aabc/parameter_particle_validation_result_set.cpp aabc/parameter_particle_validation_result_set.hpp \
				aabc/parameter_particle_validator.cpp aabc/parameter_particle_validator.hpp \
				aabc/parameter_particle_validator_multi_loci_polytable.cpp aabc/parameter_particle_validator_multi_loci_polytable.hpp \
				aabc/parameter_particle_validator_multi_loci_polytable_heterozygosity.cpp aabc/parameter_particle_validator_multi_loci_polytable_heterozygosity.hpp \
				aabc/pop_mc_abc_simulator_adaptive.cpp aabc/pop_mc_abc_simulator_adaptive.hpp \
				aabc/parameter_prior.cpp aabc/parameter_prior.hpp \
				aabc/parameter_prior_normal.cpp aabc/parameter_prior_normal.hpp \
				aabc/parameter_prior_uniform.cpp aabc/parameter_prior_uniform.hpp \
				aabc/parameter_prior_set.cpp aabc/parameter_prior_set.hpp \
				aabc/mixture_test.cpp aabc/mixture_test.hpp \
				aabc/mixture_test_sampler.cpp aabc/mixture_test_sampler.hpp \
				aabc/normalsd_test.cpp aabc/normalsd_test.hpp \
				aabc/normalsd_test_sampler.cpp aabc/normalsd_test_sampler.hpp \
				aabc/mixture_test_2d.cpp aabc/mixture_test_2d.hpp \
				aabc/mixture_test_2d_sampler.cpp aabc/mixture_test_2d_sampler.hpp \
				reg/reg_linreg_step.cpp reg/reg_linreg_step.hpp \
				reg/reg_transformer.cpp reg/reg_transformer.hpp \
				reg/reg_transformer_log.cpp reg/reg_transformer_log.hpp \
				reg/reg_transformer_tan.cpp reg/reg_transformer_tan.hpp \
				reg/reg_transformer_null.cpp reg/reg_transformer_null.hpp
				

#basic am compiler flags
AM_CPPFLAGS = $(BOOST_CPPFLAGS) @SRC_SUBS_INCLUDES@

#all the flags and make sure we have the AM_CPPFLAGS too or it will be overridden	
#libmct_a_CPPFLAGS = $(AM_CPPFLAGS)

AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) 

#LDADD = ${top_builddir}/libmscplusplus-0.0.a
LDADD = -lmscplusplus-0.0 


LIBS += $(BOOST_PROGRAM_OPTIONS_LIBS)
