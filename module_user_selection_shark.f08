! **********************************************************************************************************************************
! This module specifies selection functions of different galaxy surveys
!
! For an illustration with comments of how a selection function is coded, please look at the subroutine selection_example().
!
! To add a new selection function, follow these steps:
! 1) Copy the function selection_example() to the bottom of the code
! 2) Change the function name to the name your survey: selection_[survey name]()
! 3) Code up the selection function by completing the five 'case' clauses; use selected = .true. for clauses that are not required
! 4) Link the selection function to a survey name by adding a new 'case' clause in the subroutine assign_selection_function
!
! Notes:
! + see the example of the WALLABY survey for an illustration of multiple surveys with very similar selection functions
! **********************************************************************************************************************************


module module_user_selection

! **********************************************************************************************************************************
! MODULE INTERFACE (DO NOT EDIT)
! **********************************************************************************************************************************

use shared_module_core
use shared_module_maths
use shared_module_cosmology
use module_global
use module_parameters
use module_conversion
use module_user_routines
use module_selection_tools

private

public   :: assign_selection_function

contains

! **********************************************************************************************************************************
! LINKING SURVEY NAMES TO SELECTION FUNCTIONS
! **********************************************************************************************************************************

! All custom selection functions must be linked to a survey name, i.e. the parameter "survey" in the parameter file, via a
! case clause in the following subroutine.

subroutine assign_selection_function

   select case (trim(para%survey))
      case('example');              selection_function => selection_example
      case('gama');                 selection_function => selection_gama
      case('devils');               selection_function => selection_devils
      case('waves-wide');           selection_function => selection_waves_wide
      case('waves-deep');           selection_function => selection_waves_deep
   case default
      call selection_function_unknown
   end select   

end subroutine assign_selection_function


! **********************************************************************************************************************************
! CUSTOM SELECTION FUNCTIONS (private, only accessible via the public pointer "selection_function")
! **********************************************************************************************************************************

! Default example (do not edit this example, but use it as a template for new selection functions)

subroutine selection_example(pos,sam,sky,range,selected)

   ! do not edit
   implicit none
   type(type_spherical),intent(in),optional     :: pos      ! has components dc [length unit of parameterfile], ra [deg], dec [deg]
   type(type_sam),intent(in),optional           :: sam      ! has components as defined in the "module_user_routines_..."
   type(type_sky_galaxy),intent(in),optional    :: sky      ! has components as defined in the "module_user_routines_..."
   type(type_fov),intent(inout),optional        :: range    ! has components dc(2) [length unit], ra(2) [deg], dec(2) [deg]
   logical,intent(inout),optional               :: selected
   ! end do not edit
   
   select case (selection_type(pos,sam,sky,range,selected))
   case (return_position_range)
   
      ! here enter the individual maximal ranges of comoving distance, right ascension and declination covered by the survey,
      ! as restrictive as possible; these ranges are mandatory
      ! Note: it is possible to have the range%ra go from a larger number (e.g. 330) to a smaller one (e.g. 30). In this case,
      !       the sky defined by the wedge from 330 deg (=-30deg) to 30 deg is considered.
      range%dc = (/0.0,300.0/)      ! [simulation length units, here Mpc/h] comoving distance range
      range%ra = (/150.0,210.0/)    ! [deg] range of right ascensions, bound to 0 to 360
      range%dec = (/-10.0,10.0/)    ! [deg] range of declinations, bound to -90 to +90
      
   case (select_by_pos)
   
      ! here enter additional restrictions on the position (pos)
      ! if all positional restrictions are already covered by the ranges above, leave this clause empty
      selected = pos%ra>=190.0 .or. pos%ra<=170.0
      
   case (select_by_sam)
      
      ! here add additional, maximal restrictions that only use the SAM-properties in type type_sam;
      ! if no such restrictions exist, leave this clause empty
      selected = sam%mstars_disk+sam%mstars_bulge>1e8
      
   case (select_by_pos_and_sam)
   
      ! here add additional, maximal restrictions that require both the position (pos) *and* SAM-properties (sam);
      ! if no such restrictions exist, leave this clause empty
      selected = (sam%mstars_disk+sam%mstars_bulge)/pos%dc**2>1e4 ! rough preselection, only for acceleration
      
   case (select_by_all)
   
      ! here add additional, maximal restrictions that require apparent sky properties (sky), as defined in module_user_routines,
      ! possibly combined with position and SAM properties;
      ! if no such restrictions exist, leave this clause empty
      selected = sky%mag<19.23 .and. sky%zobs<0.1 ! select by apparent magnitude and redshift
      
   end select
   
end subroutine

! **********************************************************************************************************************************

! GAMA survey
! This function defines a pre-selection of the GAMA galaxy survey, before applying the final cut by apparent magnitude. To apply
! the final cuts the stingray output must first be post-processed by Viperfish (by A. Robotham), which generates observer-frame
! SEDs as a function of the star formation history and position of each source.

subroutine selection_gama(pos,sam,sky,range,selected)

   implicit none
   type(type_spherical),intent(in),optional     :: pos
   type(type_sam),intent(in),optional           :: sam
   type(type_sky_galaxy),intent(in),optional    :: sky
   type(type_fov),intent(inout),optional        :: range
   logical,intent(inout),optional               :: selected
   
   ! computation variables
   real*4            :: mstars ! [Msun] stellar mass
   real*4            :: dl ! [Mpc] comoving distance
   real*4            :: mag ! generic apparent magnitude assuming M/L=1
   real*4,parameter  :: dmag = 4.0 ! magnitude tolerance
   
   ! selection function
   select case (selection_type(pos,sam,sky,range,selected))
   case (return_position_range)
      range%dc = (/0.0,2450.0/)     ! [simulation length units, here Mpc/h] distance range (to z~0.6)
      range%ra = (/129.0,351.0/)    ! [deg] range of right ascensions, bound to 0 to 360
      range%dec = (/-35.0,3.0/)     ! [deg] range of declinations, bound to -90 to +90
   case (select_by_pos)
      selected = ((pos%ra> 129.0).and.(pos%ra<=141.0).and.(pos%dec>= -2.00).and.(pos%dec<= +3.00)).or. & ! field G09
               & ((pos%ra>=174.0).and.(pos%ra<=186.0).and.(pos%dec>= -3.00).and.(pos%dec<= +2.00)).or. & ! field G12
               & ((pos%ra>=211.5).and.(pos%ra<=223.5).and.(pos%dec>= -2.00).and.(pos%dec<= +3.00)).or. & ! field G15
               & ((pos%ra>=339.0).and.(pos%ra<=351.0).and.(pos%dec>=-35.00).and.(pos%dec<=-30.00))       ! field G23
   case (select_by_sam)
      selected = (sam%mstars_disk+sam%mstars_bulge)/para%h>1e7
   case (select_by_pos_and_sam)
      ! note: The selection below is based on a rough estimate of a generic apparent magnitude mag.
      !       This same magnitude is computed later and stored in sky%mag, which is used onder 'selected_by_all'.
      !       The point of pre-computing the magnitude here, using the comoving distance as an
      !       inferior limit for the luminosity distance, is to accelerate the computation by avoiding time-consuming
      !       computations of many apparent galaxy properties
      mstars = (sam%mstars_disk+sam%mstars_bulge)/para%h ! [Msun]
      dl = pos%dc/para%h ! [Mph] comoving distance as an inferior limit for the luminosity distance, which would require sky%zobs
      mag = convert_absmag2appmag(convert_stellarmass2absmag(mstars,1.0),dl)
      selected = mag<=19.65+dmag
   case (select_by_all)
      selected = sky%mag<=19.65+dmag
   end select
   
end subroutine

! **********************************************************************************************************************************

! DEVILS survey
! This function defines a pre-selection of the DEVILS galaxy survey, before applying the final cut by apparent magnitude. To apply
! the final cuts the stingray output must first be post-processed by Viperfish (by A. Robotham), which generates observer-frame
! SEDs as a function of the star formation history and position of each source.

subroutine selection_devils(pos,sam,sky,range,selected)

   implicit none
   type(type_spherical),intent(in),optional     :: pos
   type(type_sam),intent(in),optional           :: sam
   type(type_sky_galaxy),intent(in),optional    :: sky
   type(type_fov),intent(inout),optional        :: range
   logical,intent(inout),optional               :: selected
   
   ! computation variables
   real*4            :: mag ! rough estimate of a apparent magnitude assuming M/L=1
   real*4,parameter  :: dmag = 4.0 ! magnitude tolerance
   
   ! selection function
   select case (selection_type(pos,sam,sky,range,selected))
   case (return_position_range)
      range%dc = (/0.0,4000.00/) ! [simulation length units, here Mpc/h] distance range (to z~1.2)
      range%ra = (/34.0,150.70/) ! [deg] range of right ascensions, bound to 0 to 360
      range%dec = (/-28.5,2.79/) ! [deg] range of declinations, bound to -90 to +90
   case (select_by_pos)
      selected = ((pos%ra>= 34.000).and.(pos%ra<= 37.050).and.(pos%dec>= -5.200).and.(pos%dec<= -4.200)).or. & ! D02 (XMM-LSS)
               & ((pos%ra>= 52.263).and.(pos%ra<= 53.963).and.(pos%dec>=-28.500).and.(pos%dec<=-27.500)).or. & ! D03 (ECDFS)
               & ((pos%ra>=149.380).and.(pos%ra<=150.700).and.(pos%dec>= +1.650).and.(pos%dec<= +2.790))       ! D10 (COSMOS)
   case (select_by_sam)
      selected = (sam%mstars_disk+sam%mstars_bulge)/para%h>1e7
   case (select_by_pos_and_sam)
      ! note: see comments in selection_gama
      mag = convert_absmag2appmag(convert_stellarmass2absmag((sam%mstars_disk+sam%mstars_bulge)/para%h,1.0),pos%dc/para%h)
      selected = mag<=21.2+dmag
   case (select_by_all)
      selected = sky%mag<=21.2+dmag
   end select
   
end subroutine

! **********************************************************************************************************************************

! WAVES-Wide survey

subroutine selection_waves_wide(pos,sam,sky,range,selected)

   implicit none
   type(type_spherical),intent(in),optional     :: pos
   type(type_sam),intent(in),optional           :: sam
   type(type_sky_galaxy),intent(in),optional    :: sky
   type(type_fov),intent(inout),optional        :: range
   logical,intent(inout),optional               :: selected
   
   ! computation variables
   real*4            :: mag ! rough estimate of a apparent magnitude assuming M/L=1
   real*4,parameter  :: dmag = 4.0 ! magnitude tolerance
   
   ! selection function
   select case (selection_type(pos,sam,sky,range,selected))
   case (return_position_range)
      range%dc = (/0.0,1250.0/) ! [simulation length units, here Mpc/h] distance range (to z~0.3)
      range%ra = (/0.0,360.0/)  ! [deg] range of right ascensions, bound to 0 to 360
      range%dec = (/-40.0,5.0/) ! [deg] range of declinations, bound to -90 to +90
   case (select_by_pos)
      selected = ((pos%ra>= 157.000).and.(pos%ra<= 225.000).and.(pos%dec>= -3.000).and.(pos%dec<=  4.000)).or. & ! WWN
               & ((pos%ra>=   0.000).and.(pos%ra<=  52.500).and.(pos%dec>=-35.900).and.(pos%dec<=-27.000)).or. & ! WWS (East)
               & ((pos%ra>= 330.000).and.(pos%ra<  360.000).and.(pos%dec>=-35.900).and.(pos%dec<=-27.000))       ! WWS (West)
   case (select_by_sam)
      selected = (sam%mstars_disk+sam%mstars_bulge)/para%h>1e7
   case (select_by_pos_and_sam)
      ! note: see comments in selection_gama
      mag = convert_absmag2appmag(convert_stellarmass2absmag((sam%mstars_disk+sam%mstars_bulge)/para%h,1.0),pos%dc/para%h)
      selected = mag<=21.1+dmag
   case (select_by_all)
      selected = sky%mag<=21.1+dmag
   end select
   
end subroutine

! **********************************************************************************************************************************

! WAVES-Deep survey

subroutine selection_waves_deep(pos,sam,sky,range,selected)

   implicit none
   type(type_spherical),intent(in),optional     :: pos
   type(type_sam),intent(in),optional           :: sam
   type(type_sky_galaxy),intent(in),optional    :: sky
   type(type_fov),intent(inout),optional        :: range
   logical,intent(inout),optional               :: selected
   
   ! computation variables
   real*4            :: mag ! rough estimate of a apparent magnitude assuming M/L=1
   real*4,parameter  :: dmag = 4.0 ! magnitude tolerance
   
   ! selection function
   select case (selection_type(pos,sam,sky,range,selected))
   case (return_position_range)
      range%dc = (/0.0,4000.0/) ! [simulation length units, here Mpc/h] distance range (to z~1.2)
      range%ra = (/-21.0,151.12/) ! [deg] range of right ascensions, bound to 0 to 360
      range%dec = (/-44.7,3.5/) ! [deg] range of declinations, bound to -90 to +90
   case (select_by_pos)
      selected = ((pos%ra>=339.00).and.(pos%ra<=351.00).and.(pos%dec>=-35.00).and.(pos%dec<=-30.00)).or. & ! G23
               & ((pos%ra>=  7.95).and.(pos%ra<=  9.95).and.(pos%dec>=-44.70).and.(pos%dec<=-42.70)).or. & ! ELAIS-S
               & ((pos%ra>= 34.50).and.(pos%ra<= 36.50).and.(pos%dec>= -6.55).and.(pos%dec<= -4.55)).or. & ! XMM-LSS
               & ((pos%ra>= 52.00).and.(pos%ra<= 54.00).and.(pos%dec>=-29.00).and.(pos%dec<=-27.00)).or. & ! eCDFS
               & ((pos%ra>=149.12).and.(pos%ra<=151.12).and.(pos%dec>=  1.50).and.(pos%dec<=  3.50))       ! eCOSMOS
   case (select_by_sam)
      selected = (sam%mstars_disk+sam%mstars_bulge)/para%h>1e7
   case (select_by_pos_and_sam)
      ! note: see comments in selection_gama
      mag = convert_absmag2appmag(convert_stellarmass2absmag((sam%mstars_disk+sam%mstars_bulge)/para%h,1.0),pos%dc/para%h)
      selected = mag<=24.0+dmag
   case (select_by_all)
      selected = sky%mag<=24.0+dmag
   end select
   
end subroutine

! **********************************************************************************************************************************

end module module_user_selection
