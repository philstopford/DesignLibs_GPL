﻿using Burkardt.Types;

namespace Burkardt.Square
{
    public static partial class MinimalRule
    {
        public static double[] smr35()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SMR35 returns the SMR rule of degree 35.
            //
            //  Discussion:
            // 
            //    DEGREE: 35
            //    SYMMETRY: (X,  Y),  (-X, -Y)
            //    POINTS CARDINALITY: 220
            //    NORM INF MOMS. RESIDUAL: 8.88178e-16
            //    SUM NEGATIVE WEIGHTS: 0
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
            //    Output, double *SMR35[3*220], the requested rule.
            //
        {
            int degree = 35;
            int order;
            double[] xw =
            {
                -9.963054076831271e-01, 5.895227652492729e-01, 3.040448117975760e-03,
                -9.959258876623405e-01, 8.226232166208993e-01, 1.978505799274471e-03,
                -9.954067560360648e-01, -9.074145716693054e-01, 1.810898452206362e-03,
                -9.952722138287154e-01, 2.785041956973013e-01, 4.045664611880734e-03,
                -9.949775104476455e-01, -7.198973208744759e-01, 2.841977406946536e-03,
                -9.943750550352027e-01, -3.487179202856331e-02, 4.548201806973092e-03,
                -9.943677037885296e-01, 9.710410041880194e-01, 1.118097006661641e-03,
                -9.937038655576630e-01, -9.931050540649843e-01, 4.830386117047890e-04,
                -9.930202681222473e-01, -3.214092657950156e-01, 4.905950396961080e-03,
                -9.861727841676639e-01, -5.400704381496431e-01, 6.797444141494801e-03,
                -9.810071147719678e-01, 8.985653300431585e-01, 3.181435032108599e-03,
                -9.740062631571137e-01, 7.162386303873299e-01, 7.887232214456035e-03,
                -9.739587062789258e-01, 4.439901118287898e-01, 1.034505230549930e-02,
                -9.726572697503230e-01, -8.181619946621148e-01, 6.464658905526641e-03,
                -9.698067055262342e-01, -9.633285091443817e-01, 3.282158060938487e-03,
                -9.694166462805192e-01, 1.256612526452283e-01, 1.203688131505905e-02,
                -9.658674622089268e-01, 9.953452266059633e-01, 1.102717947131251e-03,
                -9.638166221951805e-01, -1.751088896720577e-01, 1.281494987287384e-02,
                -9.504683215166834e-01, 9.376314986069719e-01, 4.086139127613287e-03,
                -9.449990557642718e-01, -6.629472122415720e-01, 1.284843891927071e-02,
                -9.446169086076862e-01, -4.135452324650374e-01, 1.479906753649862e-02,
                -9.317772095678226e-01, 8.201150672577782e-01, 9.940838732765474e-03,
                -9.290486544464471e-01, -8.950190921642169e-01, 7.779964394420137e-03,
                -9.283184737295814e-01, 5.879017021514434e-01, 1.561783824168382e-02,
                -9.245462285990703e-01, 2.974342857591706e-01, 1.847140993304593e-02,
                -9.211763298375819e-01, -9.934526556341349e-01, 2.003191789183488e-03,
                -9.119151052554647e-01, -4.274927442267495e-02, 1.424779189715308e-02,
                -9.057117846078058e-01, 4.659054057320846e-02, 9.353615677783169e-03,
                -8.999935234232900e-01, 9.752903231391948e-01, 4.574852292020883e-03,
                -8.829230087690153e-01, -7.770114466508429e-01, 1.555785382570586e-02,
                -8.824467573040489e-01, -2.752146559569857e-01, 2.057326622412542e-02,
                -8.735663725309100e-01, -5.411359973847757e-01, 2.141510206682392e-02,
                -8.722180260064613e-01, 8.841005984426642e-01, 8.223856754800985e-03,
                -8.688253850828861e-01, -9.504355153693655e-01, 7.016335896084859e-03,
                -8.610711834248425e-01, 7.150978677167033e-01, 1.850564562244901e-02,
                -8.562600324120904e-01, 4.495124318073283e-01, 2.357983853078307e-02,
                -8.358095013326707e-01, 1.733093930153186e-01, 2.670539971994917e-02,
                -8.228302442461702e-01, -1.697016693626475e-01, 1.342647121783435e-02,
                -8.174385098192266e-01, 9.961711079884141e-01, 2.165794610184791e-03,
                -8.166426016882831e-01, 9.271485107406205e-01, 7.934960713327414e-03,
                -7.985152091160521e-01, -8.693626506598062e-01, 1.561754453690349e-02,
                -7.981120997224805e-01, -9.863110488956950e-01, 4.428350131300970e-03,
                -7.904587180395604e-01, -3.698232163504092e-02, 2.186616762504613e-02,
                -7.887652330718642e-01, -6.709289186517926e-01, 2.371430123275542e-02,
                -7.796834866912661e-01, -3.994715783200819e-01, 2.994104825767884e-02,
                -7.704963456494349e-01, 5.900999992322654e-01, 2.533980569560729e-02,
                -7.640326515570037e-01, 8.125887805195392e-01, 1.949174894848284e-02,
                -7.450492234413376e-01, 3.427381908986717e-01, 2.313801537268639e-02,
                -7.328136392431331e-01, 2.875348312510743e-01, 1.087089248290535e-02,
                -7.295286264509984e-01, 9.665226078227548e-01, 8.147629083103793e-03,
                -7.084311098316746e-01, -1.768747230410797e-01, 2.397513960514677e-02,
                -7.005552670374735e-01, -9.406601330339404e-01, 1.246369990823231e-02,
                -6.897531031934704e-01, 6.782883935108303e-01, 6.655878361359627e-03,
                -6.848525475182752e-01, -7.813397580356689e-01, 2.347027715648843e-02,
                -6.787362302215644e-01, 8.353105038445782e-02, 3.660462342407370e-02,
                -6.758525517932752e-01, -5.427274842488018e-01, 3.178580319617729e-02,
                -6.598249437273979e-01, -9.948642685624477e-01, 2.808772731098325e-03,
                -6.499799074317244e-01, 8.924083156022753e-01, 1.749067931565174e-02,
                -6.406772691253437e-01, 7.139228812806040e-01, 2.272673241453971e-02,
                -6.404909510129633e-01, 4.864679391303438e-01, 2.960589457868493e-02,
                -6.308221736554908e-01, -2.945267523482654e-01, 3.146868659787020e-02,
                -6.143455982549889e-01, 9.920826636566066e-01, 4.351963504883606e-03,
                -5.844710675512947e-01, 2.592109992058108e-01, 2.283706283744573e-02,
                -5.748230680261558e-01, -8.760965122066653e-01, 1.943926421467558e-02,
                -5.572474716940677e-01, -6.686420737982589e-01, 3.083669208006901e-02,
                -5.482102814487624e-01, 4.226896138717067e-01, 1.428889064276768e-02,
                -5.435737968160658e-01, -7.127399923480068e-02, 4.496806151459203e-02,
                -5.391683921956473e-01, -9.726355827179098e-01, 9.266285418087179e-03,
                -5.207225999437588e-01, 2.085864473946192e-01, 2.344644441982280e-02,
                -5.185875206855339e-01, 8.093271136372929e-01, 2.465400557427630e-02,
                -5.177821997159030e-01, -4.860941365404695e-01, 4.035705383914671e-03,
                -5.159801356223465e-01, 9.489058642290333e-01, 1.356103767266921e-02,
                -5.158924660295786e-01, -4.251904122231722e-01, 3.560324497451925e-02,
                -4.994044935070279e-01, 6.099462217677528e-01, 3.154922279377764e-02,
                -4.542210356062989e-01, -7.958468688022844e-01, 1.794582284359643e-02,
                -4.056267255080869e-01, 5.341375520421459e-01, 1.152087899889739e-02,
                -4.055948117557182e-01, -9.243438360258911e-01, 1.650133309775475e-02,
                -4.040537647708318e-01, -9.964418579498329e-01, 3.051316164369587e-03,
                -4.039882331666476e-01, -2.369652785329704e-01, 4.786645009346951e-02,
                -4.035957925466688e-01, -7.396121234399951e-01, 1.283147906403060e-02,
                -3.998139656878879e-01, 3.785894153380550e-01, 3.993632109430833e-02,
                -3.963413695053817e-01, 9.888743963929636e-01, 6.292415031278153e-03,
                -3.916192571372809e-01, 8.570203854231701e-02, 4.709957356606227e-02,
                -3.781489380730715e-01, -5.664445529989698e-01, 3.920699629566945e-02,
                -3.716587536775290e-01, 8.852446618690641e-01, 2.141062768778597e-02,
                -3.634328788865228e-01, 7.287742400150319e-01, 2.448742112847983e-02,
                -2.987755245919356e-01, 6.636636932174717e-01, 1.342406173089754e-02,
                -2.968276435525883e-01, -8.304928792595524e-01, 1.658848141268110e-02,
                -2.730446394137964e-01, -9.731854164765670e-01, 9.487914154560882e-03,
                -2.527215209939402e-01, -3.953707636595537e-01, 4.787894718529206e-02,
                -2.484030447599341e-01, 9.526421508958445e-01, 1.422151685597499e-02,
                -2.449636454507225e-01, -7.846110838738764e-02, 5.146237215684395e-02,
                -2.425768385346887e-01, 5.209707813340320e-01, 3.945757174858697e-02,
                -2.350501139206987e-01, 2.470184275136251e-01, 4.984023350983210e-02,
                -2.242013177906116e-01, -6.901076261298711e-01, 3.533647809963283e-02,
                -2.185301091710398e-01, -8.814048828689003e-01, 1.288022793213647e-02,
                -2.041253198436666e-01, 7.913052621344840e-01, 2.097344678642826e-02,
                -1.711186961614260e-01, 8.484142146354349e-01, 1.191006595550569e-02,
                -1.706560382250964e-01, 9.949901599692252e-01, 3.586815258663693e-03,
                -1.215442732714673e-01, -9.952093452521580e-01, 3.571400312739183e-03,
                -1.212739677278976e-01, -9.375172954719908e-01, 1.267403468943651e-02,
                -9.463646218871977e-02, -5.420702319082150e-01, 4.444769202139284e-02,
                -8.912243328958734e-02, 6.529263044412255e-01, 3.773327086808358e-02,
                -8.859975837369544e-02, -2.416302562834947e-01, 5.169192291225042e-02,
                -8.213238766551616e-02, -8.377296354933957e-01, 1.300368186580167e-02,
                -7.960209828358143e-02, 8.301571435730447e-02, 5.326449572668888e-02,
                -6.910083831974916e-02, 3.985144846429817e-01, 4.829901192754266e-02,
                -6.703676747477079e-02, 9.036691320042987e-01, 1.766559618844168e-02,
                -3.859244773141641e-02, -7.670969426559235e-01, 2.567105291338053e-02,
                -2.835677761623644e-02, 9.748999327824780e-01, 8.816515305029195e-03,
                2.835677761623644e-02, -9.748999327824780e-01, 8.816515305029195e-03,
                3.859244773141641e-02, 7.670969426559235e-01, 2.567105291338053e-02,
                6.703676747477079e-02, -9.036691320042987e-01, 1.766559618844168e-02,
                6.910083831974916e-02, -3.985144846429817e-01, 4.829901192754266e-02,
                7.960209828358143e-02, -8.301571435730447e-02, 5.326449572668888e-02,
                8.213238766551616e-02, 8.377296354933957e-01, 1.300368186580167e-02,
                8.859975837369544e-02, 2.416302562834947e-01, 5.169192291225042e-02,
                8.912243328958734e-02, -6.529263044412255e-01, 3.773327086808358e-02,
                9.463646218871977e-02, 5.420702319082150e-01, 4.444769202139284e-02,
                1.212739677278976e-01, 9.375172954719908e-01, 1.267403468943651e-02,
                1.215442732714673e-01, 9.952093452521580e-01, 3.571400312739183e-03,
                1.706560382250964e-01, -9.949901599692252e-01, 3.586815258663693e-03,
                1.711186961614260e-01, -8.484142146354349e-01, 1.191006595550569e-02,
                2.041253198436666e-01, -7.913052621344840e-01, 2.097344678642826e-02,
                2.185301091710398e-01, 8.814048828689003e-01, 1.288022793213647e-02,
                2.242013177906116e-01, 6.901076261298711e-01, 3.533647809963283e-02,
                2.350501139206987e-01, -2.470184275136251e-01, 4.984023350983210e-02,
                2.425768385346887e-01, -5.209707813340320e-01, 3.945757174858697e-02,
                2.449636454507225e-01, 7.846110838738764e-02, 5.146237215684395e-02,
                2.484030447599341e-01, -9.526421508958445e-01, 1.422151685597499e-02,
                2.527215209939402e-01, 3.953707636595537e-01, 4.787894718529206e-02,
                2.730446394137964e-01, 9.731854164765670e-01, 9.487914154560882e-03,
                2.968276435525883e-01, 8.304928792595524e-01, 1.658848141268110e-02,
                2.987755245919356e-01, -6.636636932174717e-01, 1.342406173089754e-02,
                3.634328788865228e-01, -7.287742400150319e-01, 2.448742112847983e-02,
                3.716587536775290e-01, -8.852446618690641e-01, 2.141062768778597e-02,
                3.781489380730715e-01, 5.664445529989698e-01, 3.920699629566945e-02,
                3.916192571372809e-01, -8.570203854231701e-02, 4.709957356606227e-02,
                3.963413695053817e-01, -9.888743963929636e-01, 6.292415031278153e-03,
                3.998139656878879e-01, -3.785894153380550e-01, 3.993632109430833e-02,
                4.035957925466688e-01, 7.396121234399951e-01, 1.283147906403060e-02,
                4.039882331666476e-01, 2.369652785329704e-01, 4.786645009346951e-02,
                4.040537647708318e-01, 9.964418579498329e-01, 3.051316164369587e-03,
                4.055948117557182e-01, 9.243438360258911e-01, 1.650133309775475e-02,
                4.056267255080869e-01, -5.341375520421459e-01, 1.152087899889739e-02,
                4.542210356062989e-01, 7.958468688022844e-01, 1.794582284359643e-02,
                4.994044935070279e-01, -6.099462217677528e-01, 3.154922279377764e-02,
                5.158924660295786e-01, 4.251904122231722e-01, 3.560324497451925e-02,
                5.159801356223465e-01, -9.489058642290333e-01, 1.356103767266921e-02,
                5.177821997159030e-01, 4.860941365404695e-01, 4.035705383914671e-03,
                5.185875206855339e-01, -8.093271136372929e-01, 2.465400557427630e-02,
                5.207225999437588e-01, -2.085864473946192e-01, 2.344644441982280e-02,
                5.391683921956473e-01, 9.726355827179098e-01, 9.266285418087179e-03,
                5.435737968160658e-01, 7.127399923480068e-02, 4.496806151459203e-02,
                5.482102814487624e-01, -4.226896138717067e-01, 1.428889064276768e-02,
                5.572474716940677e-01, 6.686420737982589e-01, 3.083669208006901e-02,
                5.748230680261558e-01, 8.760965122066653e-01, 1.943926421467558e-02,
                5.844710675512947e-01, -2.592109992058108e-01, 2.283706283744573e-02,
                6.143455982549889e-01, -9.920826636566066e-01, 4.351963504883606e-03,
                6.308221736554908e-01, 2.945267523482654e-01, 3.146868659787020e-02,
                6.404909510129633e-01, -4.864679391303438e-01, 2.960589457868493e-02,
                6.406772691253437e-01, -7.139228812806040e-01, 2.272673241453971e-02,
                6.499799074317244e-01, -8.924083156022753e-01, 1.749067931565174e-02,
                6.598249437273979e-01, 9.948642685624477e-01, 2.808772731098325e-03,
                6.758525517932752e-01, 5.427274842488018e-01, 3.178580319617729e-02,
                6.787362302215644e-01, -8.353105038445782e-02, 3.660462342407370e-02,
                6.848525475182752e-01, 7.813397580356689e-01, 2.347027715648843e-02,
                6.897531031934704e-01, -6.782883935108303e-01, 6.655878361359627e-03,
                7.005552670374735e-01, 9.406601330339404e-01, 1.246369990823231e-02,
                7.084311098316746e-01, 1.768747230410797e-01, 2.397513960514677e-02,
                7.295286264509984e-01, -9.665226078227548e-01, 8.147629083103793e-03,
                7.328136392431331e-01, -2.875348312510743e-01, 1.087089248290535e-02,
                7.450492234413376e-01, -3.427381908986717e-01, 2.313801537268639e-02,
                7.640326515570037e-01, -8.125887805195392e-01, 1.949174894848284e-02,
                7.704963456494349e-01, -5.900999992322654e-01, 2.533980569560729e-02,
                7.796834866912661e-01, 3.994715783200819e-01, 2.994104825767884e-02,
                7.887652330718642e-01, 6.709289186517926e-01, 2.371430123275542e-02,
                7.904587180395604e-01, 3.698232163504092e-02, 2.186616762504613e-02,
                7.981120997224805e-01, 9.863110488956950e-01, 4.428350131300970e-03,
                7.985152091160521e-01, 8.693626506598062e-01, 1.561754453690349e-02,
                8.166426016882831e-01, -9.271485107406205e-01, 7.934960713327414e-03,
                8.174385098192266e-01, -9.961711079884141e-01, 2.165794610184791e-03,
                8.228302442461702e-01, 1.697016693626475e-01, 1.342647121783435e-02,
                8.358095013326707e-01, -1.733093930153186e-01, 2.670539971994917e-02,
                8.562600324120904e-01, -4.495124318073283e-01, 2.357983853078307e-02,
                8.610711834248425e-01, -7.150978677167033e-01, 1.850564562244901e-02,
                8.688253850828861e-01, 9.504355153693655e-01, 7.016335896084859e-03,
                8.722180260064613e-01, -8.841005984426642e-01, 8.223856754800985e-03,
                8.735663725309100e-01, 5.411359973847757e-01, 2.141510206682392e-02,
                8.824467573040489e-01, 2.752146559569857e-01, 2.057326622412542e-02,
                8.829230087690153e-01, 7.770114466508429e-01, 1.555785382570586e-02,
                8.999935234232900e-01, -9.752903231391948e-01, 4.574852292020883e-03,
                9.057117846078058e-01, -4.659054057320846e-02, 9.353615677783169e-03,
                9.119151052554647e-01, 4.274927442267495e-02, 1.424779189715308e-02,
                9.211763298375819e-01, 9.934526556341349e-01, 2.003191789183488e-03,
                9.245462285990703e-01, -2.974342857591706e-01, 1.847140993304593e-02,
                9.283184737295814e-01, -5.879017021514434e-01, 1.561783824168382e-02,
                9.290486544464471e-01, 8.950190921642169e-01, 7.779964394420137e-03,
                9.317772095678226e-01, -8.201150672577782e-01, 9.940838732765474e-03,
                9.446169086076862e-01, 4.135452324650374e-01, 1.479906753649862e-02,
                9.449990557642718e-01, 6.629472122415720e-01, 1.284843891927071e-02,
                9.504683215166834e-01, -9.376314986069719e-01, 4.086139127613287e-03,
                9.638166221951805e-01, 1.751088896720577e-01, 1.281494987287384e-02,
                9.658674622089268e-01, -9.953452266059633e-01, 1.102717947131251e-03,
                9.694166462805192e-01, -1.256612526452283e-01, 1.203688131505905e-02,
                9.698067055262342e-01, 9.633285091443817e-01, 3.282158060938487e-03,
                9.726572697503230e-01, 8.181619946621148e-01, 6.464658905526641e-03,
                9.739587062789258e-01, -4.439901118287898e-01, 1.034505230549930e-02,
                9.740062631571137e-01, -7.162386303873299e-01, 7.887232214456035e-03,
                9.810071147719678e-01, -8.985653300431585e-01, 3.181435032108599e-03,
                9.861727841676639e-01, 5.400704381496431e-01, 6.797444141494801e-03,
                9.930202681222473e-01, 3.214092657950156e-01, 4.905950396961080e-03,
                9.937038655576630e-01, 9.931050540649843e-01, 4.830386117047890e-04,
                9.943677037885296e-01, -9.710410041880194e-01, 1.118097006661641e-03,
                9.943750550352027e-01, 3.487179202856331e-02, 4.548201806973092e-03,
                9.949775104476455e-01, 7.198973208744759e-01, 2.841977406946536e-03,
                9.952722138287154e-01, -2.785041956973013e-01, 4.045664611880734e-03,
                9.954067560360648e-01, 9.074145716693054e-01, 1.810898452206362e-03,
                9.959258876623405e-01, -8.226232166208993e-01, 1.978505799274471e-03,
                9.963054076831271e-01, -5.895227652492729e-01, 3.040448117975760e-03
            };
            double[] xw_copy;

            order = square_minimal_rule_order(degree);
            xw_copy = typeMethods.r8mat_copy_new(3, order, xw);

            return xw_copy;
        }
    }
}