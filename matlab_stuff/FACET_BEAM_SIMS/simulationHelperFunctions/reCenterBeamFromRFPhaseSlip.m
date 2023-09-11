function beam = reCenterBeamFromRFPhaseSlip(beam)
% Track to end of each bending section and re-center beam in each case to
% take care of phase-slip with respect to RF and orbit excursions due to SR
% energy losses. Also remove linear dispersion induced by CSR effects.
% In reality this is done with RF phasing, orbit feedbacks and beam tuning.
  for ii=1:5
    if ii<5
      coef=polyfit(beam.Bunch.x(6,~beam.Bunch.stop),beam.Bunch.x(ii,~beam.Bunch.stop),1);
      coef(2)=0;
      beam.Bunch.x(ii,~beam.Bunch.stop)=beam.Bunch.x(ii,~beam.Bunch.stop)-polyval(coef,beam.Bunch.x(6,~beam.Bunch.stop));
    end
    beam.Bunch.x(ii,~beam.Bunch.stop)=beam.Bunch.x(ii,~beam.Bunch.stop)-mean(beam.Bunch.x(ii,~beam.Bunch.stop));
  end