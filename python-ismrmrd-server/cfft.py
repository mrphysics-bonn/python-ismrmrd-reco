
try:
    import pyfftw.interfaces.numpy_fft as fft
except ImportError:
    import numpy.fft as fft


def cfftn(data, axes):
    """ Centered fast fourier transform, n-dimensional.

    :param data: Complex input data.
    :param axes: Axes along which to shift and transform.
    :return: Fourier transformed data.
    """
    return fft.fftshift(fft.fftn(fft.ifftshift(data, axes=axes), axes=axes, norm='ortho'), axes=axes)


def cifftn(data, axes):
    """ Centered inverse fast fourier transform, n-dimensional.

    :param data: Complex input data.
    :param axes: Axes along which to shift.
    :return: Inverse fourier transformed data.
    """
    return fft.ifftshift(fft.ifftn(fft.fftshift(data, axes=axes), axes=axes, norm='ortho'), axes=axes)
