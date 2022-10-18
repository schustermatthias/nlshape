import numpy as np

quadrature_rules = {
    '16points': [np.array([[0.33333333, 0.33333333],
                           [0.45929259, 0.45929259],
                           [0.45929259, 0.08141482],
                           [0.08141482, 0.45929259],
                           [0.17056931, 0.17056931],
                           [0.17056931, 0.65886138],
                           [0.65886138, 0.17056931],
                           [0.05054723, 0.05054723],
                           [0.05054723, 0.89890554],
                           [0.89890554, 0.05054723],
                           [0.26311283, 0.72849239],
                           [0.72849239, 0.00839478],
                           [0.00839478, 0.26311283],
                           [0.72849239, 0.26311283],
                           [0.26311283, 0.00839478],
                           [0.00839478, 0.72849239]]),
                 np.array([0.14431560767779,
                           0.09509163426728,
                           0.09509163426728,
                           0.09509163426728,
                           0.10321737053472,
                           0.10321737053472,
                           0.10321737053472,
                           0.03245849762320,
                           0.03245849762320,
                           0.03245849762320,
                           0.02723031417443,
                           0.02723031417443,
                           0.02723031417443,
                           0.02723031417443,
                           0.02723031417443,
                           0.02723031417443])
                 ],
    '3points': [np.array([[1. / 6., 1. / 6.],
                          [1. / 6., 2. / 3.],
                          [2. / 3., 1. / 6.]]),
                1. / 3 * np.ones(3)
                ],
    '4points': [np.array([[1. / 3., 1. / 3.],
                          [0.2, 0.6],
                          [0.2, 0.2],
                          [0.6, 0.2]]),
                np.array([-27. / 48.,
                          25. / 48.,
                          25. / 48.,
                          25. / 48.])
                ]
}
Px = np.array([[0.33333333333333, 0.33333333333333],
               [0.47014206410511, 0.47014206410511],
               [0.47014206410511, 0.05971587178977],
               [0.05971587178977, 0.47014206410511],
               [0.10128650732346, 0.10128650732346],
               [0.10128650732346, 0.79742698535309],
               [0.79742698535309, 0.10128650732346]])
dx = 0.5 * np.array([0.22500000000000,
                     0.13239415278851,
                     0.13239415278851,
                     0.13239415278851,
                     0.12593918054483,
                     0.12593918054483,
                     0.12593918054483])
Py = Px
dy = dx
configuration = {
    'nu': 0.001,

    'save_results': 1,

    'l_bfgs': 1,
    'memory_length': 5,
    'armijo_const': 0.0001,

    'target_shape': 'target_shape',
    'init_shape': 'example_1',

    'number_iterations': 20,

    'c_per': 0.1,
    'remesh': 1,
    'element_size': 0.05,

    'boundary_label': 3,
    'interface_label': 12,

    'alpha': 0.5,

    'adapt_mu_max': 0,
    'lambda': 0.0,
    'mu_min': 0.0,
    'mu_max': 0.5,
    'lame_threshold_up': 4,
    'lame_threshold_down': 1,
    'adapt_up': 1.2,
    'adapt_down': 0.8,

    'source': [100, -10],

    'kernel_dependent': 1,
    #fractional_shape, integrable_unsym
    'kernel': dict(function="integrable_unsym", horizon=1. / 10., outputdim=1, fractional_s=0.5),
    'nlfem_conf': dict(ansatz="CG",
                       approxBalls={"method": "retriangulate",
                                    "isPlacePointOnCap": True,  # required for "retriangulate" only
                                    # "averageBallWeights": [1., 1., 1.]  # required for "averageBall" only
                                    },
                       closeElements="fractional",
                       quadrature={"outer": {"points": Px,
                                             "weights": dx
                                             },
                                   "inner": {"points": Py,
                                             "weights": dy
                                             },
                                   "tensorGaussDegree": 5  # Degree of tensor Gauss quadrature for singular kernels.
                                   },
                       is_ShapeDerivative=0,
                       verbose=False),
    'nlfem_shape_conf': dict(ansatz="CG",  # only CG possible
                          approxBalls={"method": "retriangulate_shape",
                                       "isPlacePointOnCap": True,  # required for "retriangulate" only
                                        # "averageBallWeights": [1., 1., 1.]  # required for "averageBall" only
                                       },
                          closeElements="fractional_shape",
                          quadrature={"outer": {"points": Px,
                                                "weights": dx
                                                },
                                      "inner": {"points": Py,
                                                "weights": dy
                                                },
                                      "tensorGaussDegree": 5  # Degree of tensor Gauss quadrature for singular kernels.
                                      },
                          is_ShapeDerivative=1,
                          verbose=False)
}
# END
