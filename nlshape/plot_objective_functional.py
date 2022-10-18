import matplotlib.pyplot as plt

ex1_sym = [0.026196677427335474, 0.02408817646855967, 0.021299913780789125, 0.017858425704825476, 0.01397926910273219,
           0.01014500226287106, 0.0071473270539797115, 0.0022596051859988713, 0.001557808812353154,
           0.001535381374382983, 0.0016836124406080703, 0.0016283547446950299, 0.0015754434751357879,
           0.0015743974720475345, 0.001571234949296336, 0.001569264211830409, 0.0015685503249294182,
           0.001567958293973684, 0.0015679435047172752, 0.0015679336077894796, 0.0015679268570491366]
ex2_sym = [0.019956318894254133, 0.01695483922870268, 0.013727799677679965, 0.01055432374067951, 0.007552132229067797,
           0.006067759683788461, 0.004509016710587695, 0.001987006598368243, 0.0017588352674291555,
           0.0016923057002552779, 0.0017460015468286642, 0.0017187988880177035, 0.0016484773721480357,
           0.0016304943445995333, 0.001594985203400116, 0.0015773563351075737, 0.0015754936752985014,
           0.0015749063959135889, 0.0015749049610518093, 0.001574902885909556, 0.0015749004226671478]
ex1_non = [0.07261763086867132, 0.05367371727232279, 0.023561678188234708, 0.020239311690505198, 0.011270414364587793,
           0.007375893167831072, 0.001770533248947149, 0.0015847154639358684, 0.001580158534617124,
           0.0015795397655350935, 0.0015711099974534402, 0.0015704797376827535, 0.0015703913273177631,
           0.001570011389317983, 0.0015692309453340773, 0.0015681235150843587, 0.001567777850064216,
           0.0015672003471950337, 0.0015670669515354008, 0.0015670361713009861, 0.0015670283226765676]
ex2_non = [0.0555787507114254, 0.04047629968861391, 0.032194258310373335, 0.028238184577266588, 0.01966021799295683,
           0.016517742682784536, 0.0028060030026334846, 0.0018669012128756381, 0.001719123486228345,
           0.0017097758960876533, 0.0016817174595114943, 0.0016759973461433408, 0.0015902568567052423,
           0.0015861877654057505, 0.0015845813930331137, 0.0015798718161807751, 0.001575308993729038,
           0.0015679649749005037, 0.0015674358686694674, 0.0015671390548553544, 0.0015671248172839261]

plt.plot(ex1_sym, label='Example 1 singular symmetric kernel')
plt.plot(ex2_sym, label='Example 2 singular symmetric kernel')
plt.plot(ex1_non, label='Example 1 square integrable kernel')
plt.plot(ex2_non, label='Example 2 square integrable kernel')

plt.xlim([0, 20])
plt.xlabel('number of iterations')
plt.ylabel('objective functional value')
plt.legend()
plt.savefig('comparison_objective_functional.pdf')
plt.show()
