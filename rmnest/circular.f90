

subroutine circular(frequency,&
                    frequency_cen,&
                    len_ch,&
                    psi_0,&
                    grm,&
                    alpha,&
                    chi,&
                    phi,&
                    theta,&
                    rot_vector)
        
        implicit none

        !f2py threadsafe
        integer(kind=4), intent(in):: len_ch
                
        integer(kind=4):: i,j

        !Defining constants
        
        real(kind=8):: const_c
        real(kind=8):: const_mega
        real(kind=8):: deg2rad
        

        real(kind=8), intent(in):: frequency_cen
        real(kind=8), intent(in):: alpha
        real(kind=8), intent(in):: chi
        real(kind=8), intent(in):: phi
        real(kind=8), intent(in):: grm
        real(kind=8), intent(in):: psi_0
        real(kind=8), intent(in):: theta
        
        real(kind=8), intent(in), dimension(len_ch)       ::frequency
        real(kind=8),  dimension(len_ch)       ::psi

        real(kind=8),  dimension(len_ch, 3)       ::stokesp
        real(kind=8),  dimension(3, 3)            ::rot_kernel
        real(kind=8), intent(out), dimension(len_ch, 3)         ::rot_vector
        
        const_c     =   299792458
        const_mega  =   1000000

        do i=1, len_ch
                psi(i)  =   psi_0*deg2rad + grm * (&
                             ((const_c/(frequency(i) * const_mega))**alpha)&
                             -((const_c/(frequency_cen * const_mega))**alpha))

                stokesp(i,1)   =   cos(2*psi(i))*cos(2*chi*deg2rad)
                stokesp(i,2)   =   sin(2*psi(i))*cos(2*chi*deg2rad)
                stokesp(i,3)   =   sin(2*chi*deg2rad) 
        end do
        
        !do i=1, 3 
        !    do j=1, len_ch
        !        sq(i)(j)   =   cos(2*psi(i))*cos(2*chi*deg2rad)
        !        su(i)(j)   =   sin(2*psi(i))*cos(2*chi*deg2rad)
        !        sv(i)(j)   =   sin(2*chi*deg2rad) 
        !    end do
        !end do
        
        !defining manualy :')
        rot_kernel(1,1)    =   cos(theta*deg2rad)*cos(phi*deg2rad)
        rot_kernel(1,2)    =  -cos(theta*deg2rad)*sin(phi*deg2rad)
        rot_kernel(1,3)    =   sin(theta*deg2rad)
        rot_kernel(2,1)    =   sin(phi*deg2rad)
        rot_kernel(2,2)    =   cos(phi*deg2rad)
        rot_kernel(2,3)    =      0 
        rot_kernel(3,1)    =   -sin(theta*deg2rad)*cos(phi*deg2rad)
        rot_kernel(3,2)    =   sin(theta*deg2rad)*sin(phi*deg2rad)
        rot_kernel(3,3)    =   cos(theta*deg2rad)
        
        !mat_thea = array([[cos(theta)*cos(phi), -cos(theta)*sin(phi), sin(theta)],
        ![sin(phi), cos(phi), 0], [-sin(theta)*cos(phi), sin(theta)*sin(phi),
        !cos(theta)]]) 
        
        rot_vector          =   matmul(stokesp, rot_kernel)

        end subroutine circular

