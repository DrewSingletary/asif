#ifndef _ASIF_LEARNING_UTILS_H_
#define _ASIF_LEARNING_UTILS_H_

#include "asif_utils.h"

namespace ASIF
{
	typedef struct
	{
		uint32_t d_drift_in = 0;
		uint32_t d_act_in = 0;
		uint32_t d_drift_hidden = 0;
		uint32_t d_act_hidden = 0;
		uint32_t d_drift_hidden_2 = 0;
		uint32_t d_act_hidden_2 = 0;
		uint32_t d_drift_out = 0;
		uint32_t d_act_out = 0;
		const double* w_1_drift;
		const double* w_2_drift;
		const double* w_3_drift;
		const double* b_1_drift;
		const double* b_2_drift;
		const double* b_3_drift;
		const double* w_1_act;
		const double* w_2_act;
		const double* w_3_act;
		const double* b_1_act;
		const double* b_2_act;
		const double* b_3_act;
		double* Lfh_diff;
		double* Lgh_diff;
	} LearningData;

	inline void driftNN(LearningData *data_, const double *drift_input, double *drift_output)
		{
		    double output_1[data_->d_drift_hidden]; // Define output of first NN layer
		    // Multiply input by first weight matrix
		    matrixVectorMultiply(data_->w_1_drift, data_->d_drift_hidden, data_->d_drift_in,
		                         drift_input, data_->d_drift_in,
		                         output_1);

		    // Add bias and use ReLu nonlinearity for each output of first weight matrix
		    for(uint32_t i=0;i<data_->d_drift_hidden;i++)
		    {
		      // Add bias, then use max with 0 to compute ReLu
		      output_1[i] = fmax(0., output_1[i] + data_->b_1_drift[i]);
		    }

		    double output_2[data_->d_drift_hidden_2]; // Define output of first NN layer
				// Multiply input by second weight matrix
				matrixVectorMultiply(data_->w_2_drift, data_->d_drift_hidden_2, data_->d_drift_hidden,
														 output_1, data_->d_drift_hidden,
														 output_2);

				// Add bias and use ReLu nonlinearity for each output of first weight matrix
				for(uint32_t i=0;i<data_->d_drift_hidden_2;i++)
				{
					// Add bias, then use max with 0 to compute ReLu
					output_2[i] = fmax(0., output_2[i] + data_->b_2_drift[i]);
				}

		    double output_3[data_->d_drift_out]; // Define output of second NN layer

		    // Multiply output of first layer by second weight matrix.
		    matrixVectorMultiply(data_->w_3_drift, data_->d_drift_out, data_->d_drift_hidden_2,
		                         output_2, data_->d_drift_hidden_2,
		                         output_3);

		    // Add bias for each output of the second weight matrix
		    for(uint32_t i=0;i<data_->d_drift_out;i++)
		    {
		      // Add bias
		      drift_output[i] = output_3[i] + data_->b_3_drift[i];
		    }
		}

	inline void actNN(LearningData *data_, const double *act_input, double *act_output)
		{
		double output_1[data_->d_act_hidden]; // Define output of first NN layer
				    // Multiply input by first weight matrix
				    matrixVectorMultiply(data_->w_1_act, data_->d_act_hidden, data_->d_act_in,
				                         act_input, data_->d_act_in,
				                         output_1);

				    // Add bias and use ReLu nonlinearity for each output of first weight matrix
				    for(uint32_t i=0;i<data_->d_act_hidden;i++)
				    {
				      // Add bias, then use max with 0 to compute ReLu
				      output_1[i] = fmax(0., output_1[i] + data_->b_1_act[i]);
				    }

				    double output_2[data_->d_act_hidden_2]; // Define output of first NN layer
						// Multiply input by second weight matrix
						matrixVectorMultiply(data_->w_2_act, data_->d_act_hidden_2, data_->d_act_hidden,
																 output_1, data_->d_act_hidden,
																 output_2);

						// Add bias and use ReLu nonlinearity for each output of first weight matrix
						for(uint32_t i=0;i<data_->d_act_hidden_2;i++)
						{
							// Add bias, then use max with 0 to compute ReLu
							output_2[i] = fmax(0., output_2[i] + data_->b_2_act[i]);
						}

				    double output_3[data_->d_act_out]; // Define output of second NN layer

				    // Multiply output of first layer by second weight matrix.
				    matrixVectorMultiply(data_->w_3_act, data_->d_act_out, data_->d_act_hidden_2,
				                         output_2, data_->d_act_hidden_2,
				                         output_3);

				    // Add bias for each output of the second weight matrix
				    for(uint32_t i=0;i<data_->d_act_out;i++)
				    {
				      // Add bias
				      act_output[i] = output_3[i] + data_->b_3_act[i];
				    }

	}

	inline void update_weights(LearningData *data_, const double *x, const uint32_t nx, const double *Dh, double *Lfh, double *Lgh, const uint32_t nu) {
		double drift_input[data_->d_drift_in] = {0.};
		for (int i = 0; i < (int) nx; i ++) {
			drift_input[i] = x[i];
		}

		for (int i = 0; i < (int) nx; i++) {
			drift_input[i+nx] = Dh[i];
		}

		double act_input[data_->d_act_in] = {0.};
		for (int i = 0; i < (int) nx; i ++) {
			act_input[i] = x[i];
		}
		for (int i = 0; i < (int) nx; i++) {
			act_input[i+nx] = Dh[i];
		}

		double drift_output[data_->d_drift_out] = {0.};
		double act_output[data_->d_act_out] = {0.};

		driftNN(data_, drift_input, drift_output);
		actNN(data_, act_input, act_output);

		data_->Lfh_diff = new double(drift_output[0]);
		data_->Lgh_diff = new double[nu];

		Lfh[0]+=drift_output[0];

		for(uint32_t i=0;i<nu;i++)
		{
			Lgh[i]+=act_output[i];
			data_->Lgh_diff[i] = act_output[i];
		}
	}

} //end ASIF namespace

#endif
