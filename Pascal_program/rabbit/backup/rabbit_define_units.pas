unit rabbit_define_units;

{$mode ObjFPC}{$H+}

interface

uses
  Classes, SysUtils;

type           {here you declare the data structure for you individuals}
  Array2DInt = array of array of integer;
  Array3DInt = array of array of array of integer;
  MapOfLists = array of array of Tlist;

  PAgent = ^Agent;

  Agent = record
    sex: string[1];
    Age: integer;      // months
    Pregnant: boolean;  // Individual pregnant?
  end;

  var

    RabbitPopulationSpatial: MapOfLists;
    RabbitPopulation: TList;
    Rabbit: PAgent;
    n_ini: integer;
    rep_prob: real;
    AlphaR: real;
    BetaR: real;
    surv_prob: real;
    AlphaS: real;
    BetaS: real;
    sett_prob: real;
    avg_steps: integer;
    AlphaD: real;
    BetaD: real;
    sink: real;
    RabbitPopulationSize, RabbitPopulationSpatialSize: longint;
    n_sim: integer;
    n_extint: integer;
    sum_distance_X: integer;
    sum_distance_Y: integer;
    current_sim: integer;

    BreedingMonths_xy: Array3DInt;
    DryMonths_xy: Array3DInt;

    mapname_rabbit, breedname, dryname: string;
    dx: array[0..8] of integer = (0, 0, 1, 1, 1, 0, -1, -1, -1);
    dy: array[0..8] of integer = (0, 1, 1, 0, -1, -1, -1, 0, 1);
    tempX, tempY: integer;
    nbadmonths: byte;
    max_years: integer;
    day: integer;
    Kdensity, KCapacity, kCapacity_low, kCapacity_high: real;
    lambda, dens_opt, sigma: real;                            // Dispersal parameters
    max_age, juv_age, sub_adult_age, adult_age: integer;      // Ages in months
    disp_age: integer;
    MeanLitSize, SdLitSize: real;
    r_int, r_dens_effect, r_second_effect, r_later_effect: real;
    MortP_at_month_old, s_int, s_extra_juv,s_food_effect, s_dens_effect: real;
    check_disp: array of integer;
    check_disp_count: integer;
    mcmc_file_location: string;

  const
    //diames: array[0..11] of integer = (1, 31, 62, 92, 122, 153, 183, 214, 244, 274, 305, 335);
    //diarep: array[0..11] of integer = (15, 45, 75, 105, 135, 165, 195, 225, 255, 285, 315, 345);
    days_in_month: array[1..12] of integer = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);

    disp_x: array[0..48] of integer = (-3, -2, -1, 0, 1, 2, 3,
                                       -3, -2, -1, 0, 1, 2, 3,
                                       -3, -2, -1, 0, 1, 2, 3,
                                       -3, -2, -1, 0, 1, 2, 3,
                                       -3, -2, -1, 0, 1, 2, 3,
                                       -3, -2, -1, 0, 1, 2, 3,
                                       -3, -2, -1, 0, 1, 2, 3);
    disp_y: array[0..48] of integer = (3, 3, 3, 3, 3, 3, 3,
                                       2, 2, 2, 2, 2, 2, 2,
                                       1, 1, 1, 1, 1, 1, 1,
                                       0, 0, 0, 0, 0, 0, 0,
                                      -1, -1, -1, -1, -1, -1, -1,
                                      -2, -2, -2, -2, -2, -2, -2,
                                      -3, -3, -3, -3, -3, -3, -3);
    disp_dist: array[0..48] of integer = (3, 3, 3, 3, 3, 3, 3,
                                          3, 2, 2, 2, 2, 2, 3,
                                          3, 2, 1, 1, 1, 2, 3,
                                          3, 2, 1, 0, 1, 2, 3,
                                          3, 2, 1, 1, 1, 2, 3,
                                          3, 2, 2, 2, 2, 2, 3,
                                          3, 3, 3, 3, 3, 3, 3);

    MCMC_PARAMETER_SELECTION_RUN: boolean = false;

implementation

end.

