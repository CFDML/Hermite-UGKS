program test_temperature_convention
use GASTHEORY, only: get_lambda_from_temperature, get_temperature_from_lambda
implicit none

call require(abs(get_lambda_from_temperature(1.1d0) - 1.d0/1.1d0) < 1.d-15, &
             'lambda is not the inverse dimensionless temperature')
call require(abs(get_temperature_from_lambda(0.5d0) - 2.d0) < 1.d-15, &
             'temperature is not the inverse of lambda')

write(*,*) 'PASS: temperature and lambda use the paper convention'

contains

subroutine require(condition, message)
  logical, intent(in) :: condition
  character(len=*), intent(in) :: message
  if (.not. condition) then
    write(*,*) 'FAIL: ', trim(message)
    error stop 1
  end if
end subroutine require

end program test_temperature_convention
