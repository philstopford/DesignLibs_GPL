﻿using Burkardt.Types;

namespace Burkardt.Square
{
    public static partial class MinimalRule
    {
        public static double[] smr40()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SMR40 returns the SMR rule of degree 40.
            //
            //  Discussion:
            // 
            //    DEGREE: 40
            //    POINTS CARDINALITY: 296
            //    NORM INF MOMS. RESIDUAL: 8.88178e-16
            //    SUM NEGATIVE WEIGHTS: 0.00000e+00,  
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 February 2018
            //
            //  Author:
            //
            //    Original MATLAB version by Mattia Festa, Alvise Sommariva,
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Mattia Festa, Alvise Sommariva,
            //    Computing almost minimal formulas on the square,
            //    Journal of Computational and Applied Mathematics,
            //    Volume 17, Number 236, November 2012, pages 4296-4302.
            //
            //  Parameters:
            //
            //    Output, double *SMR40[3*296], the requested rule.
            //
        {
            int degree = 40;
            int order;
            double[] xw =
            {
                -9.980037227196614e-01, 1.255939443038006e-01, 2.053102393463664e-03,
                -9.975374042816840e-01, 5.936108405943805e-01, 1.276187489810420e-03,
                -9.972861466402049e-01, -6.627439420820114e-01, 1.580298156253981e-03,
                -9.968592673741854e-01, 8.736423351683437e-01, 1.202945592410809e-03,
                -9.966589451189311e-01, -4.369849380251707e-01, 2.115586380569520e-03,
                -9.965219690837533e-01, -1.944638806631471e-01, 2.429919318140257e-03,
                -9.964358393204891e-01, 4.092031891829934e-01, 2.182777026867640e-03,
                -9.963941450446951e-01, 9.775574537230374e-01, 5.781357112514488e-04,
                -9.955864247426816e-01, -9.784933334191470e-01, 6.423485659734820e-04,
                -9.955242411687316e-01, -8.201261852647251e-01, 1.452148584658119e-03,
                -9.913130830294992e-01, 7.451950336474148e-01, 3.099740554076827e-03,
                -9.900337734476057e-01, -9.112486663047842e-01, 2.026525375041507e-03,
                -9.844939696556546e-01, -3.396948505159310e-02, 6.297322525864257e-03,
                -9.828708241253864e-01, 2.672014458471039e-01, 6.110796337925116e-03,
                -9.813352997502048e-01, 9.317282602187734e-01, 2.388338374639091e-03,
                -9.812506993249306e-01, -5.599844336596295e-01, 5.433547359777527e-03,
                -9.811994633761572e-01, 5.954700474743184e-01, 3.855922634832198e-03,
                -9.787064347309005e-01, -3.245284003636001e-01, 6.507740619311661e-03,
                -9.780647700362322e-01, -7.464363616852314e-01, 4.386694918546395e-03,
                -9.768709333873531e-01, 9.964868257322258e-01, 5.743979761358160e-04,
                -9.743127159036291e-01, -9.969203399240063e-01, 5.885958513598154e-04,
                -9.719890138043383e-01, 4.635529385440017e-01, 5.543426475209034e-03,
                -9.663961744212400e-01, 8.297161334267321e-01, 5.377245930028534e-03,
                -9.624486391679733e-01, -9.552403379650487e-01, 2.935271926976710e-03,
                -9.571101155486035e-01, -8.504016909299451e-01, 5.548045896834978e-03,
                -9.569508072526900e-01, 1.145259922791730e-01, 1.019364727373913e-02,
                -9.539620291105575e-01, 9.656118875216714e-01, 2.000317967473577e-03,
                -9.531361804534030e-01, -1.756012259335100e-01, 1.012788272053494e-02,
                -9.503256987566595e-01, 6.828228592655408e-01, 8.213806957862843e-03,
                -9.455690174869477e-01, -4.553679978362609e-01, 1.002546929027743e-02,
                -9.444080894387320e-01, -6.600309149596145e-01, 8.342165274488325e-03,
                -9.363352068572794e-01, 3.305989409491193e-01, 1.063637549554008e-02,
                -9.344111368509483e-01, 5.018191785468655e-01, 7.093229113292136e-03,
                -9.250445674046065e-01, 8.992422617422033e-01, 6.279305248680532e-03,
                -9.225703808148057e-01, 9.858292479280879e-01, 1.834573979737579e-03,
                -9.165441617582062e-01, -9.853128908447532e-01, 2.565324347525059e-03,
                -9.150614770382599e-01, -3.139530918232870e-02, 1.359896871759167e-02,
                -9.103915600727824e-01, -9.092647756707973e-01, 6.216724896190497e-03,
                -9.074945980525453e-01, -3.110157725349993e-01, 1.298683840693059e-02,
                -9.043535181463653e-01, -7.713160602214696e-01, 1.024303332784389e-02,
                -9.025325558225301e-01, 7.743255817267981e-01, 1.011918005042377e-02,
                -8.936138114540279e-01, 1.960899239677014e-01, 1.411969841594799e-02,
                -8.921757562307675e-01, 5.847775228961032e-01, 1.223494929054215e-02,
                -8.914413155620955e-01, -5.644405075920796e-01, 1.287350206631372e-02,
                -8.651606789959365e-01, 9.507227302652309e-01, 5.817570712124294e-03,
                -8.624964059071404e-01, 4.007509803672114e-01, 1.631963389227952e-02,
                -8.624016491888394e-01, 9.973607974975222e-01, 1.085292256544525e-03,
                -8.571481676731489e-01, -1.732871980674100e-01, 1.643819257648864e-02,
                -8.511908836822574e-01, -4.127628928140689e-01, 1.185185065054727e-02,
                -8.471275338153255e-01, -9.533530423203328e-01, 5.827596397572629e-03,
                -8.433676706517138e-01, 5.983417524756374e-02, 1.569370700632547e-02,
                -8.424649352939378e-01, 8.531606755621376e-01, 1.000254175836821e-02,
                -8.382814166059266e-01, -9.973371948408554e-01, 1.236406010086881e-03,
                -8.380565845515963e-01, -8.466246757209853e-01, 1.074234765384687e-02,
                -8.321936110008199e-01, -6.779831075876902e-01, 1.522781703701599e-02,
                -8.286278033800711e-01, 6.833732733978219e-01, 1.509304227664999e-02,
                -8.055255758765038e-01, -4.763880300923605e-01, 1.085366933632746e-02,
                -7.971673647784905e-01, 2.532323209821991e-01, 2.052836615250460e-02,
                -7.910441928184301e-01, 5.077971678500216e-01, 1.780077542770354e-02,
                -7.859172990177037e-01, -7.833283965104983e-02, 1.472595970793188e-02,
                -7.798224683613482e-01, 9.817158588687893e-01, 4.411372664168903e-03,
                -7.733447858053083e-01, -2.912095449947633e-01, 2.109987317572403e-02,
                -7.707288564403509e-01, 9.135259150232076e-01, 8.345329225642115e-03,
                -7.578198426154767e-01, -9.796098169168449e-01, 4.456574433684158e-03,
                -7.562414896721606e-01, -9.049214717524879e-01, 1.001882084525143e-02,
                -7.548211869038036e-01, 7.785365721991255e-01, 1.478039679943451e-02,
                -7.494378413286128e-01, -7.662608330090768e-01, 1.581956505815484e-02,
                -7.445473618861123e-01, -5.699700714382581e-01, 1.889632922544271e-02,
                -7.423038822674627e-01, 7.897961507665216e-02, 1.715604114405427e-02,
                -7.170452369551987e-01, 6.041689391836269e-01, 1.721434133773814e-02,
                -7.096177258303671e-01, 3.669023638814550e-01, 2.303538633973196e-02,
                -6.947187172157958e-01, -9.098243292500578e-02, 1.600840192685054e-02,
                -6.827474344539016e-01, 8.632122010101693e-01, 9.488912002746175e-03,
                -6.766766916048634e-01, -4.118520712505950e-01, 2.365726668658341e-02,
                -6.757537207484740e-01, 9.965428761249577e-01, 1.913682375287724e-03,
                -6.716977448629430e-01, 9.507129459344548e-01, 7.775380155384150e-03,
                -6.705672963967729e-01, -2.084540920008020e-01, 1.628326394109379e-02,
                -6.634191394897307e-01, -9.956567399175887e-01, 2.049561637770933e-03,
                -6.628050205606234e-01, 1.672916166402376e-01, 2.265410339187529e-02,
                -6.572284494011310e-01, -9.451593820143880e-01, 8.033722247688203e-03,
                -6.502802198668451e-01, -6.663478504843525e-01, 2.086896333372951e-02,
                -6.500944652911421e-01, -8.388046887319990e-01, 1.533892815116460e-02,
                -6.494395416059597e-01, 6.978896106401927e-01, 1.632586135380591e-02,
                -6.182525423851862e-01, 4.772666066191493e-01, 2.314379901273175e-02,
                -6.097521926338272e-01, 8.178387432803899e-01, 1.027570595259021e-02,
                -5.857038932505656e-01, -1.750697876397254e-03, 2.826026275637352e-02,
                -5.765612790071245e-01, -5.184987079049800e-01, 2.276145801578116e-02,
                -5.737142051029732e-01, -3.132765352046541e-01, 1.443899895963975e-02,
                -5.644865970292597e-01, 9.783893306945404e-01, 5.710986844724012e-03,
                -5.644553463342945e-01, 2.782921931969836e-01, 2.781479815231922e-02,
                -5.580769546100564e-01, -9.759275455445522e-01, 5.684058667060024e-03,
                -5.542741288051631e-01, 9.016977074485676e-01, 1.168589067769203e-02,
                -5.542073588689838e-01, -2.043031348070334e-01, 1.827083660577868e-02,
                -5.433087051306172e-01, 5.922649560095551e-01, 1.990937282981727e-02,
                -5.386567969465166e-01, -8.983031172510674e-01, 1.272427456995568e-02,
                -5.373382383966947e-01, -7.520144164271236e-01, 2.091503656277224e-02,
                -5.224304686125945e-01, 7.455130741439847e-01, 1.698684695502466e-02,
                -4.893827154330137e-01, -5.941538126001544e-01, 1.598966080652984e-02,
                -4.722069619254669e-01, -3.587550371799506e-01, 2.132063146049633e-02,
                -4.718832732594200e-01, -9.958106507002302e-01, 2.226108918584357e-03,
                -4.706670553696944e-01, 1.154506383897263e-01, 2.935766164731613e-02,
                -4.626454569472719e-01, 9.963250210519803e-01, 2.223046715737653e-03,
                -4.600264076712823e-01, 3.966582665227552e-01, 2.629091387627286e-02,
                -4.481481045273694e-01, 9.395445338170469e-01, 8.839978177696708e-03,
                -4.415878841168007e-01, -1.062065683181053e-01, 3.058245417110486e-02,
                -4.412583635452328e-01, 5.087399017780769e-01, 1.125492379082715e-02,
                -4.304912794659402e-01, -4.469668566457143e-01, 1.012062172404055e-02,
                -4.287211566152106e-01, 8.321571485937996e-01, 1.786605784125267e-02,
                -4.223793646176035e-01, -8.333753189124466e-01, 1.730213472188374e-02,
                -4.222588890220841e-01, -9.452641667145504e-01, 1.007216794776104e-02,
                -4.083597732182007e-01, 6.509919047599022e-01, 2.181242775531822e-02,
                -3.998205203687450e-01, -6.594190487482490e-01, 2.190933309032257e-02,
                -3.695693944824219e-01, 2.193101954125042e-01, 2.327576101611239e-02,
                -3.676144523897173e-01, 9.682350478689863e-01, 4.815232646045845e-03,
                -3.591049460581930e-01, -4.614598488900264e-01, 1.852799716750301e-02,
                -3.447647402371574e-01, -2.381537435610105e-01, 2.808527790039856e-02,
                -3.341217099662899e-01, -9.826197273812773e-01, 5.494849259173287e-03,
                -3.093195936880301e-01, 3.287238905146597e-01, 2.185303504037471e-02,
                -3.072048398998126e-01, 9.859602190406442e-01, 3.831362623751753e-03,
                -3.058002206640513e-01, 2.103076331772531e-02, 3.459647562394031e-02,
                -3.044867199555938e-01, 5.186120611581733e-01, 2.833979791519306e-02,
                -3.028681432998989e-01, -7.600469172213820e-01, 1.964686881716537e-02,
                -3.024605574449728e-01, 8.931659138204667e-01, 1.538434812191329e-02,
                -3.011532601990519e-01, 7.437794313731125e-01, 2.273787326209627e-02,
                -2.887440306182171e-01, -8.957912044678301e-01, 1.555524233275056e-02,
                -2.644411240015258e-01, -5.502788852260921e-01, 2.671171545840637e-02,
                -2.400609781144449e-01, -3.260570838043731e-01, 2.569267849167105e-02,
                -2.246489878361267e-01, -9.978033236073471e-01, 1.770873804828288e-03,
                -2.066616281756470e-01, -1.045383048279050e-01, 2.497609050660571e-02,
                -2.009374136079489e-01, 2.040282941161825e-01, 2.267428085593906e-02,
                -1.909926939589495e-01, 9.968668571166585e-01, 2.126801518422638e-03,
                -1.887172076094760e-01, -9.549438196754605e-01, 1.035110388694925e-02,
                -1.885639548532559e-01, -6.847997847346543e-01, 2.014138051533820e-02,
                -1.856269263888624e-01, 9.475270956402012e-01, 1.146925671119005e-02,
                -1.790760401308231e-01, 6.312120892250483e-01, 2.910873684757638e-02,
                -1.707661353723164e-01, 8.213933035442971e-01, 2.057181928201535e-02,
                -1.702074855482027e-01, 4.002072179476484e-01, 2.982548255683739e-02,
                -1.505878349304818e-01, -8.293254383680554e-01, 2.076026265668495e-02,
                -1.380074232134180e-01, 1.170336739663703e-01, 2.527234775085595e-02,
                -1.268389913099443e-01, -4.236318257465341e-01, 3.117645114412499e-02,
                -1.153288220374788e-01, -1.747226323716413e-01, 2.609156073209279e-02,
                -8.202442049825323e-02, -6.034256289197292e-01, 2.293044125103543e-02,
                -6.778022692220509e-02, -9.864359427377531e-01, 5.537492701577899e-03,
                -5.134750275993507e-02, 9.806699464341012e-01, 7.025195533489713e-03,
                -4.796731115368809e-02, 5.109816289271156e-01, 2.949281042166952e-02,
                -4.605994163927313e-02, 8.940716194171991e-01, 1.696854743218261e-02,
                -4.447210172626447e-02, -9.096553765331403e-01, 1.548395021034487e-02,
                -4.115113842285135e-02, 7.297743920654538e-01, 2.604976478474774e-02,
                -2.744836558183677e-02, 4.417586914136311e-03, 3.387227419183413e-02,
                -1.995488830453244e-02, 2.898374142768337e-01, 3.548548223572283e-02,
                -9.393841375026718e-03, -7.464248084281863e-01, 2.510074937831712e-02,
                5.604522124288533e-03, -2.838153304903738e-01, 3.557931287068862e-02,
                3.379282258563019e-02, -5.108284849354834e-01, 2.741093593275501e-02,
                7.396931495071729e-02, 6.103364495915078e-01, 2.553255759952799e-02,
                7.811073060419052e-02, -9.596104105748841e-01, 9.711338256956472e-03,
                8.053272615641659e-02, 1.616603676436590e-01, 2.424025867338976e-02,
                8.332876758512328e-02, -9.980434804933414e-01, 1.698291109393467e-03,
                8.611110621478242e-02, 9.972454574738652e-01, 2.305464957382294e-03,
                9.098018561909829e-02, 9.458495772991818e-01, 1.224095474595298e-02,
                9.177090470118981e-02, 8.205555321043025e-01, 2.227240136214740e-02,
                9.706026488340756e-02, -8.463187115614078e-01, 2.019012230112588e-02,
                9.815256847947500e-02, -9.239456937873544e-02, 2.644315338527209e-02,
                1.222602596423681e-01, 4.213511759613903e-01, 3.341374966024566e-02,
                1.323057369500668e-01, -6.509117254364495e-01, 2.871043880147142e-02,
                1.476438086930435e-01, -4.017405050343198e-01, 2.837624979369559e-02,
                1.700850694113563e-01, -1.780332539520100e-01, 2.291185695253523e-02,
                1.736805569574461e-01, 1.096612012416960e-01, 2.405540809236990e-02,
                1.843432890535553e-01, 6.990228911926680e-01, 2.163148537165614e-02,
                1.995648926700597e-01, -9.818731842028169e-01, 3.370075215701449e-03,
                2.081596610552691e-01, 2.925076496248445e-01, 1.903382067011115e-02,
                2.211797307297267e-01, -9.164694403878136e-01, 1.468906732639756e-02,
                2.277074094032160e-01, 9.805932530165899e-01, 7.084751090497078e-03,
                2.292021628595871e-01, 8.912694976681395e-01, 1.721353195733831e-02,
                2.353608895201955e-01, -3.491298395101170e-01, 9.892666536128697e-03,
                2.385004221656652e-01, -7.689128083096237e-01, 2.393352337395165e-02,
                2.613881314139023e-01, 5.471551577512984e-01, 3.040265099653250e-02,
                2.662205097670882e-01, -9.919548226778696e-01, 2.707867154217599e-03,
                2.680565017395204e-01, -5.400402346502616e-01, 3.088207638947784e-02,
                2.833408004935626e-01, -6.718319274714699e-03, 3.511070135380980e-02,
                2.921989695479881e-01, 7.730713136429667e-01, 1.950015697008031e-02,
                2.932051140793899e-01, -2.568476496223008e-01, 2.600581679576460e-02,
                3.010766919534110e-01, 3.879280619333716e-01, 1.093142774025551e-02,
                3.103865035997045e-01, 2.499428455332039e-01, 2.499586360508493e-02,
                3.585933886911649e-01, -8.576820377675525e-01, 1.842882832290949e-02,
                3.596570955415763e-01, -9.615962408411522e-01, 9.068153972939736e-03,
                3.618929539106814e-01, 9.449875484526220e-01, 1.169739415667232e-02,
                3.672563754612045e-01, 9.968578675032447e-01, 2.310464434967522e-03,
                3.689758663333720e-01, -6.723055595123385e-01, 2.553001008106885e-02,
                3.871364360677614e-01, -3.963308613870815e-01, 2.492421215216349e-02,
                3.982469712949650e-01, 4.487336172829021e-01, 2.309608167942591e-02,
                3.983155249621386e-01, 6.561427847972009e-01, 2.730674172786235e-02,
                4.034252589810011e-01, 8.423386731412055e-01, 1.707487237256047e-02,
                4.102541856930609e-01, -1.410736134292944e-01, 3.151312043028231e-02,
                4.168183088396271e-01, 1.263343252035515e-01, 3.270130716222084e-02,
                4.356104174313109e-01, -9.957058788304909e-01, 2.363557214215902e-03,
                4.405318740860549e-01, -4.926398721055529e-01, 1.608280706132078e-02,
                4.649771622601932e-01, 3.483599503333080e-01, 1.925449126504078e-02,
                4.709325171484977e-01, -7.695942184085469e-01, 1.705050367502311e-02,
                4.850638966079184e-01, -9.216775010693871e-01, 1.252775237505845e-02,
                4.932957374397012e-01, 9.795648836635775e-01, 6.463044419371160e-03,
                5.128665631816022e-01, 9.033981066102003e-01, 1.292621422635387e-02,
                5.150863190443300e-01, -5.934248074499449e-01, 2.363213862387783e-02,
                5.156373648667533e-01, 5.651053767086885e-01, 1.423781492860490e-02,
                5.181912316625248e-01, -2.922194304405071e-01, 2.989079134667543e-02,
                5.251717070662635e-01, 7.552962710143433e-01, 2.242429087194216e-02,
                5.269679036105911e-01, -2.177852596986115e-02, 2.944698631679846e-02,
                5.359327347133740e-01, -8.200294790454760e-01, 8.860113200433273e-03,
                5.388079818555839e-01, -9.812499958704087e-01, 4.680745828233062e-03,
                5.434833238763666e-01, 2.380283686509251e-01, 2.395538344613452e-02,
                5.496379512343986e-01, 5.253805284097529e-01, 1.632157824790180e-02,
                6.051493656005358e-01, -4.503578461535495e-01, 2.612344895845866e-02,
                6.124277898034318e-01, 9.967024781187570e-01, 2.002070456520778e-03,
                6.140603542011213e-01, 3.676905093247132e-01, 1.890620847606754e-02,
                6.149248157512081e-01, -8.737165849180063e-01, 1.296667270713843e-02,
                6.154648019220623e-01, -7.039292064085437e-01, 2.124810385786090e-02,
                6.216520374448096e-01, 9.489344195184967e-01, 8.530620092099012e-03,
                6.262996333546856e-01, -9.603966937474196e-01, 6.014411655091729e-03,
                6.266507568536149e-01, -1.753112127213321e-01, 2.574430846684179e-02,
                6.320748250590934e-01, 9.236778723276985e-02, 2.509516970454762e-02,
                6.371103207708532e-01, 8.421673270141240e-01, 1.637617879409531e-02,
                6.445697045015427e-01, 6.589930041457702e-01, 2.309339250212301e-02,
                6.611641817999792e-01, -9.971250985273443e-01, 1.780949246753485e-03,
                6.827225489763870e-01, 4.582089303489358e-01, 1.967115820258447e-02,
                7.007607922685678e-01, -5.793001008608837e-01, 2.187824067436787e-02,
                7.074131809621955e-01, -3.372902789247548e-01, 2.360006012075187e-02,
                7.117739087420435e-01, 2.241359340497962e-01, 2.258162720143967e-02,
                7.151943805770620e-01, -7.271602060694378e-02, 2.063320130758209e-02,
                7.155680922709143e-01, -7.994435082716885e-01, 1.620969075703751e-02,
                7.187947908611272e-01, -9.289105584767267e-01, 8.674348771692817e-03,
                7.188580147804288e-01, 9.801870566125678e-01, 4.556848768699447e-03,
                7.407850806335647e-01, 9.092174747735176e-01, 1.069960175297740e-02,
                7.414535252324259e-01, 7.627992252808989e-01, 1.713716583993931e-02,
                7.620781689289542e-01, 5.682929391541718e-01, 1.992594154989019e-02,
                7.635939132299792e-01, -9.831481015790724e-01, 4.153559323873574e-03,
                7.797638735837149e-01, 4.665676142002305e-02, 1.954090589438721e-02,
                7.867620038671331e-01, 3.403489048139350e-01, 1.972374555702004e-02,
                7.883553505223821e-01, -6.893416107380125e-01, 1.663667575800614e-02,
                7.942764870013759e-01, -4.707153685105653e-01, 1.983755770035392e-02,
                7.952669529963022e-01, -2.273770055235256e-01, 2.028000921211288e-02,
                7.990872541216331e-01, 9.962955158284175e-01, 1.408927698874655e-03,
                8.066178878802995e-01, -8.762511423208362e-01, 1.090696275036013e-02,
                8.274679999131711e-01, 8.480202028380042e-01, 1.137803861452648e-02,
                8.297408153069072e-01, 9.581886918512165e-01, 5.948632265683609e-03,
                8.377680145222012e-01, 6.811297101045124e-01, 1.536534914781161e-02,
                8.427716005452982e-01, -9.571708221478029e-01, 4.952439619855199e-03,
                8.488271945024087e-01, 1.666337693793623e-01, 1.887388216602816e-02,
                8.525662030669592e-01, 4.592708043695782e-01, 1.653498201879920e-02,
                8.561315417874497e-01, -9.598113535151943e-02, 1.680778035994472e-02,
                8.606876078808835e-01, -9.969922689109865e-01, 1.288211684632280e-03,
                8.613581542207239e-01, -7.784610468370953e-01, 1.096644364207657e-02,
                8.690885595988450e-01, -5.894733063271762e-01, 1.503337345041933e-02,
                8.723930444088900e-01, -3.621107058578862e-01, 1.641608892157977e-02,
                8.962187117350728e-01, 9.900514429185334e-01, 2.177703335659953e-03,
                8.978901908775218e-01, 9.150513717072788e-01, 6.528697758436882e-03,
                8.995842578268963e-01, -9.234484518613306e-01, 5.137110006320958e-03,
                9.026936084906036e-01, 7.808928092755504e-01, 1.012526399216284e-02,
                9.070441620900885e-01, 3.009963311241410e-01, 1.501588480324425e-02,
                9.077878820329093e-01, 5.794983066365614e-01, 1.246940196306337e-02,
                9.111551364648940e-01, -8.412973077435701e-01, 5.641113054429737e-03,
                9.142518109635371e-01, 2.169517484726380e-02, 1.439386845586924e-02,
                9.173217644911268e-01, -2.254212884222193e-01, 1.167670360384378e-02,
                9.282555961436508e-01, -6.912520433990101e-01, 9.933776930685268e-03,
                9.284877628036031e-01, -9.817037571629690e-01, 2.694281259701618e-03,
                9.334452240614641e-01, -4.859399349804389e-01, 1.184281101417954e-02,
                9.468439730471343e-01, 9.657765158699502e-01, 3.003223815105551e-03,
                9.506625037743694e-01, 4.376710347741515e-01, 1.025469033794625e-02,
                9.520896121093332e-01, 8.643710868609885e-01, 5.555972363422187e-03,
                9.523271992044501e-01, 6.913007415265950e-01, 7.990950689265556e-03,
                9.545272982243453e-01, -8.811333722531620e-01, 3.828862340862868e-03,
                9.554958506512122e-01, 1.588668533845889e-01, 1.058860881322785e-02,
                9.598844841990601e-01, -3.304786335896162e-01, 6.696028105510322e-03,
                9.612092716134895e-01, -1.252637834195337e-01, 8.918983273981967e-03,
                9.664864057745631e-01, -7.771023857020345e-01, 5.108888381040824e-03,
                9.696658884223082e-01, -9.486617819335378e-01, 2.629252771609383e-03,
                9.714461152005505e-01, 9.966817132006599e-01, 6.382221125358519e-04,
                9.756072026722676e-01, -9.972566545264508e-01, 5.380433266276414e-04,
                9.756270675128569e-01, -5.974850205967676e-01, 6.836267790083512e-03,
                9.794527164689676e-01, 5.643324681844055e-01, 5.410131191715718e-03,
                9.826866520761600e-01, 7.911127382228441e-01, 3.928246838882010e-03,
                9.831392239644022e-01, 9.303396888012979e-01, 2.423589086603377e-03,
                9.832274416154096e-01, 3.033661616270960e-01, 6.310915040526618e-03,
                9.843450993838505e-01, -4.150221010386896e-01, 4.095407382162835e-03,
                9.849932890167565e-01, 7.572227591204869e-03, 5.802629795760786e-03,
                9.909963827975042e-01, -2.639507972140410e-01, 3.395140293721511e-03,
                9.912765561948801e-01, -8.332891526969999e-01, 1.953926311905468e-03,
                9.934327301490941e-01, 6.231882587178860e-01, 9.287993595683975e-04,
                9.937448936136122e-01, -9.092710501130503e-01, 1.256142182103206e-03,
                9.948882230658729e-01, -7.050697257712721e-01, 2.287960809703529e-03,
                9.951704470679726e-01, 9.801745622873680e-01, 6.475410324288955e-04,
                9.952340236934687e-01, -9.783971899965007e-01, 6.447054298301821e-04,
                9.963792308208721e-01, 6.981957876294371e-01, 1.512427084884079e-03,
                9.970264968246741e-01, 4.457495638028096e-01, 2.107128717380615e-03,
                9.977348846785328e-01, -1.427747693558082e-01, 1.552913600795280e-03,
                9.977388598284106e-01, 1.545000376629158e-01, 2.038498731858248e-03,
                9.982891584585788e-01, 8.710523206892243e-01, 9.529077858582162e-04,
                9.991656973291163e-01, -5.021799411346501e-01, 1.444993876112607e-03
            };
            double[] xw_copy;

            order = square_minimal_rule_order(degree);
            xw_copy = typeMethods.r8mat_copy_new(3, order, xw);

            return xw_copy;
        }
    }
}