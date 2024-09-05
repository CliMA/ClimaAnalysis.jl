using Test
import ClimaAnalysis: Utils
import Dates

@testset "Regexp" begin

    @test Utils.match_nc_filename("bob") === nothing

    @test Utils.match_nc_filename("ta_1d_average.nc") ==
          Tuple(["ta", "1d", "average"])

    @test Utils.match_nc_filename("ta_1m_40s_inst.nc") ==
          Tuple(["ta", "1m_40s", "inst"])

    @test Utils.match_nc_filename("pfull_6.0min_max.nc") ==
          Tuple(["pfull", "6.0min", "max"])

    @test Utils.match_nc_filename("hu_inst.nc") ==
          Tuple(["hu", nothing, "inst"])
end

@testset "Squeeze" begin
    @test Utils.squeeze([[1 2] [3 4]]) == [1, 2, 3, 4]
end

@testset "nearest_index" begin
    @test Utils.nearest_index([10, 20, 30], 0) == 1
    @test Utils.nearest_index([10, 20, 30], 100) == 3
    @test Utils.nearest_index([10, 20, 30], 20) == 2
    @test Utils.nearest_index([10, 20, 30], typemax(Int)) == 3
    @test Utils.nearest_index([10, 20, 30], -typemax(Int)) == 1
end

@testset "kwargs" begin
    @test values(Utils.kwargs(a = 1, b = 2)) == (; a = 1, b = 2)
end

@testset "seconds_to_prettystr" begin
    @test Utils.seconds_to_prettystr(1) == "1s"
    @test Utils.seconds_to_prettystr(2) == "2s"
    @test Utils.seconds_to_prettystr(60) == "1m"
    @test Utils.seconds_to_prettystr(61) == "1m 1s"
    @test Utils.seconds_to_prettystr(3600) == "1h"
    @test Utils.seconds_to_prettystr(3666) == "1h 1m 6s"
    @test Utils.seconds_to_prettystr(86400) == "1d"
    @test Utils.seconds_to_prettystr(86401) == "1d 1s"
    @test Utils.seconds_to_prettystr(86522) == "1d 2m 2s"
    @test Utils.seconds_to_prettystr(24 * 60 * 60 * 365) == "1y"
    @test Utils.seconds_to_prettystr(24 * 60 * 60 * 365 + 60) == "1y 1m"
    @test Utils.seconds_to_prettystr(24 * 60 * 60 * 365 + 68) == "1y 1m 8s"

    @test Utils.seconds_to_prettystr(24 * 60 * 60 * 365 + 68.5) == "1y 1m 8.5s"
end

@testset "format_title" begin
    @test Utils.warp_string("") == ""
    @test Utils.warp_string("test", max_width = 4) == "test"
    @test Utils.warp_string("     test   ", max_width = 4) == "test"
    @test Utils.warp_string("   test") == "test"
    @test Utils.warp_string("test1", max_width = 4) == "test1"
    @test Utils.warp_string("     test1   ", max_width = 4) == "test1"
    @test Utils.warp_string("test blah", max_width = 4) == "test\nblah"
    @test Utils.warp_string("test1 test2 test3", max_width = 4) ==
          "test1\ntest2\ntest3"
    @test Utils.warp_string("abc def", max_width = 3) == "abc\ndef"
    @test Utils.warp_string("is a test", max_width = 4) == "is a\ntest"
    @test Utils.warp_string("a b c d", max_width = 2) == "a\nb\nc\nd"
    @test Utils.warp_string("a b c d e f", max_width = 5) == "a b c\nd e f"
    @test Utils.warp_string("a\tb\nc\vd\fe\rf", max_width = 11) == "a b c d e f"
end

@testset "split by season" begin
    empty_dates = Vector{Dates.DateTime}()
    @test Utils.split_by_season(empty_dates) == ([], [], [], [])

    date = [Dates.DateTime(2015, 4, 13)]
    @test Utils.split_by_season(date) ==
          ([Dates.DateTime(2015, 4, 13)], [], [], [])

    dates = [
        Dates.DateTime(2015, 1, 13),
        Dates.DateTime(2018, 2, 13),
        Dates.DateTime(1981, 7, 6),
        Dates.DateTime(1993, 11, 19),
        Dates.DateTime(2040, 4, 1),
        Dates.DateTime(2000, 8, 18),
    ]

    expected_dates = (
        [Dates.DateTime(2040, 4, 1)],
        [Dates.DateTime(1981, 7, 6), Dates.DateTime(2000, 8, 18)],
        [Dates.DateTime(1993, 11, 19)],
        [Dates.DateTime(2015, 1, 13), Dates.DateTime(2018, 2, 13)],
    )

    @test Utils.split_by_season(dates) == expected_dates
end

@testset "equispaced" begin
    equispaced = Utils._isequispaced([1.0, 2.0, 3.0])
    @test equispaced == true

    equispaced = Utils._isequispaced([0.0, 2.0, 3.0])
    @test equispaced == false
end
